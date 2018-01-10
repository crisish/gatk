package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AnnotatedVariantProducer;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvDiscoveryDataBundle;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvDiscoveryUtils;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvType;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AssemblyContigWithFineTunedAlignments;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFHeaderLines;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVVCFWriter;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;


public final class InsDelVariantDetector implements VariantDetectorFromLocalAssemblyContigAlignments {

    @Override
    public void inferSvAndWriteVCF(final JavaRDD<AssemblyContigWithFineTunedAlignments> assemblyContigs,
                                   final SvDiscoveryDataBundle svDiscoveryDataBundle){

        final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast = svDiscoveryDataBundle.referenceSequenceDictionaryBroadcast;
        final String outputPath = svDiscoveryDataBundle.outputPath;
        final Logger toolLogger = svDiscoveryDataBundle.toolLogger;

        // TODO: 11/23/17 take insertion mappings from the input and add them to VC
        final JavaPairRDD<byte[], List<ChimericAlignment>> contigSeqAndChimeras =
                assemblyContigs
                        .map( AssemblyContigWithFineTunedAlignments::getSourceContig )
                        .mapToPair(tig -> convertAlignmentIntervalToChimericAlignment(tig,
                                StructuralVariationDiscoveryArgumentCollection.
                                        DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD,
                                StructuralVariationDiscoveryArgumentCollection.
                                        DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.DEFAULT_MIN_ALIGNMENT_LENGTH,
                                referenceSequenceDictionaryBroadcast.getValue()));

        final List<VariantContext> annotatedInsDels = produceVariantsFromSimpleChimeras(contigSeqAndChimeras, svDiscoveryDataBundle);

        SVVCFWriter.writeVCF(annotatedInsDels, outputPath,
                referenceSequenceDictionaryBroadcast.getValue(), toolLogger);
    }

    /**
     * Very similar to {@link ChimericAlignment#parseOneContig(AlignedContig, SAMSequenceDictionary, boolean, int, int, boolean)}, except that
     * badly mapped (MQ < 60) 1st alignment is no longer skipped.
     */
    private static Tuple2<byte[], List<ChimericAlignment>> convertAlignmentIntervalToChimericAlignment (final AlignedContig contig,
                                                                                                        final int mapQualThresholdInclusive,
                                                                                                        final int minAlignmentBlockSize,
                                                                                                        final SAMSequenceDictionary referenceDictionary) {

        final List<AlignmentInterval> alignmentIntervals = contig.alignmentIntervals;
        final Iterator<AlignmentInterval> iterator = alignmentIntervals.iterator();
        AlignmentInterval current = iterator.next();
        final List<ChimericAlignment> results = new ArrayList<>(alignmentIntervals.size() - 1);
        final List<String> insertionMappings = new ArrayList<>();
        while ( iterator.hasNext() ) {
            final AlignmentInterval next = iterator.next();
            if (ChimericAlignment.nextAlignmentMayBeInsertion(current, next, mapQualThresholdInclusive, minAlignmentBlockSize, true)) {
                if (iterator.hasNext()) {
                    insertionMappings.add(next.toPackedString());
                    continue;
                } else {
                    break;
                }
            }

            final ChimericAlignment ca = new ChimericAlignment(current, next, insertionMappings, contig.contigName, referenceDictionary);
            final boolean validStateForThisPath = ca.isNeitherSimpleTranslocationNorIncompletePicture();
            if ( ! validStateForThisPath )
                throw new GATKException.ShouldNeverReachHereException("Mapped assembled contigs are sent down the wrong path: " +
                        "contig suggesting \"translocation\" is sent down the insert/deletion path.\n" + contig.toString());
            results.add(ca);
        }
        return new Tuple2<>(contig.contigSequence,results);
    }

    public static List<VariantContext> produceVariantsFromSimpleChimeras(final JavaPairRDD<byte[], List<ChimericAlignment>> contigSeqAndChimeras,
                                                                         final SvDiscoveryDataBundle svDiscoveryDataBundle) {

        final Broadcast<ReferenceMultiSource> referenceBroadcast = svDiscoveryDataBundle.referenceBroadcast;
        final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast = svDiscoveryDataBundle.referenceSequenceDictionaryBroadcast;
        final List<SVInterval> assembledIntervals = svDiscoveryDataBundle.assembledIntervals;
        final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast = svDiscoveryDataBundle.cnvCallsBroadcast;
        final StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection discoverStageArgs = svDiscoveryDataBundle.discoverStageArgs;
        final String sampleId = svDiscoveryDataBundle.sampleId;
        final Logger toolLogger = svDiscoveryDataBundle.toolLogger;

        final JavaPairRDD<NovelAdjacencyReferenceLocations, Iterable<ChimericAlignment>> narlsAndSources =
                contigSeqAndChimeras
                        .flatMapToPair(tigSeqAndChimeras ->
                                discoverNovelAdjacencyFromChimericAlignments(tigSeqAndChimeras,
                                        referenceSequenceDictionaryBroadcast.getValue()))   // a filter-passing contig's alignments may or may not produce novel adjacency, hence flatmap
                        .groupByKey();                                                      // group the same novel adjacency produced by different contigs together

        narlsAndSources.cache();

        SvDiscoveryUtils.evaluateIntervalsAndNarls(assembledIntervals, narlsAndSources, referenceSequenceDictionaryBroadcast,
                discoverStageArgs, toolLogger);

        List<VariantContext> annotatedVariants =
                inferTypeAndAnnotateByAlignmentSignatures(narlsAndSources, cnvCallsBroadcast, sampleId,
                        referenceBroadcast, referenceSequenceDictionaryBroadcast);

        narlsAndSources.unpersist();

        return annotatedVariants;
    }

    //==================================================================================================================

    /**
     * Given contig alignments, emit novel adjacency not present on the reference to which the locally-assembled contigs were aligned.
     */
    private static Iterator<Tuple2<NovelAdjacencyReferenceLocations, ChimericAlignment>>
    discoverNovelAdjacencyFromChimericAlignments(final Tuple2<byte[], List<ChimericAlignment>> tigSeqAndChimeras, final SAMSequenceDictionary referenceDictionary) {
        return Utils.stream(tigSeqAndChimeras._2)
                .map(ca -> new Tuple2<>(new NovelAdjacencyReferenceLocations(ca, tigSeqAndChimeras._1, referenceDictionary), ca)).iterator();
    }


    // TODO: 1/10/18 to be deprecated
    private static List<VariantContext> inferTypeAndAnnotateByAlignmentSignatures(final JavaPairRDD<NovelAdjacencyReferenceLocations, Iterable<ChimericAlignment>> narlsAndSources,
                                                                                  final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast,
                                                                                  final String sampleId,
                                                                                  final Broadcast<ReferenceMultiSource> referenceBroadcast,
                                                                                  final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast) {

        return narlsAndSources
                .mapToPair(noveltyAndEvidence -> inferType(noveltyAndEvidence._1, noveltyAndEvidence._2))                     // type inference based on novel adjacency and evidence alignments
                .map(noveltyTypeAndEvidence ->
                        annotateVariant(                                                                                      // annotate the novel adjacency and inferred type
                                noveltyTypeAndEvidence._1,
                                noveltyTypeAndEvidence._2._1,
                                null,
                                noveltyTypeAndEvidence._2._2,
                                referenceBroadcast,
                                referenceSequenceDictionaryBroadcast,
                                cnvCallsBroadcast,
                                sampleId))
                .collect();
    }

    /**
     * Given input novel adjacency and evidence chimeric alignments, infer type of variant.
     */
    private static Tuple2<NovelAdjacencyReferenceLocations, Tuple2<SvType, Iterable<ChimericAlignment>>> inferType(
            final NovelAdjacencyReferenceLocations novelAdjacency,
            final Iterable<ChimericAlignment> chimericAlignments) {
        return new Tuple2<>(novelAdjacency,
                new Tuple2<>(SvTypeInference.inferFromNovelAdjacency(novelAdjacency), chimericAlignments));
    }

    /**
     * Produces annotated variant as described in {@link GATKSVVCFHeaderLines}, given input arguments.
     */
    static VariantContext annotateVariant(final NovelAdjacencyReferenceLocations novelAdjacency,
                                          final SvType inferredType,
                                          final byte[] altHaplotypeSeq,
                                          final Iterable<ChimericAlignment> chimericAlignments,
                                          final Broadcast<ReferenceMultiSource> broadcastReference,
                                          final Broadcast<SAMSequenceDictionary> broadcastSequenceDictionary,
                                          final Broadcast<SVIntervalTree<VariantContext>> broadcastCNVCalls,
                                          final String sampleId)
            throws IOException {
        return AnnotatedVariantProducer
                .produceAnnotatedVcFromInferredTypeAndRefLocations(novelAdjacency.leftJustifiedLeftRefLoc,
                        novelAdjacency.leftJustifiedRightRefLoc.getStart(), novelAdjacency.complication,
                        inferredType, altHaplotypeSeq, chimericAlignments,
                        broadcastReference, broadcastSequenceDictionary, broadcastCNVCalls, sampleId);
    }

}
