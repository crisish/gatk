package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.EvidenceTargetLink;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.ReadMetadata;
import org.broadinstitute.hellbender.tools.spark.sv.utils.PairedStrandedIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVIntervalTree;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.List;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection;


public final class SvDiscoveryDataBundle {

    public final String sampleId;
    public String outputPath;
    public final DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection discoverStageArgs;

    public final JavaRDD<GATKRead> assemblyRawAlignments;

    public final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast;
    public final List<SVInterval> assembledIntervals;
    public final PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTargetLinks;
    public final ReadMetadata metadata;

    public final Broadcast<SAMFileHeader> headerBroadcast;
    public final Broadcast<ReferenceMultiSource> referenceBroadcast;
    public final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast;

    public final Logger toolLogger;

    public SvDiscoveryDataBundle(final JavaSparkContext ctx,
                                 final DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection discoverStageArgs,
                                 final String outputPath,
                                 final ReadMetadata metadata,
                                 final List<SVInterval> assembledIntervals,
                                 final PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTargetLinks,
                                 final JavaRDD<GATKRead> reads,
                                 final SAMFileHeader headerForReads,
                                 final ReferenceMultiSource reference,
                                 final Logger toolLogger) {

        this(SVUtils.getSampleId(headerForReads), discoverStageArgs, outputPath, metadata, assembledIntervals,
                evidenceTargetLinks, reads, toolLogger,
                ctx.broadcast(reference), ctx.broadcast(headerForReads.getSequenceDictionary()), ctx.broadcast(headerForReads),
                broadcastCNVCalls(ctx, headerForReads, SVUtils.getSampleId(headerForReads), discoverStageArgs));
    }

    public SvDiscoveryDataBundle(final String sampleId,
                                 final DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection discoverStageArgs,
                                 final String outputPath,
                                 final ReadMetadata metadata,
                                 final List<SVInterval> assembledIntervals,
                                 final PairedStrandedIntervalTree<EvidenceTargetLink> evidenceTargetLinks,
                                 final JavaRDD<GATKRead> reads,
                                 final Logger toolLogger,
                                 final Broadcast<ReferenceMultiSource> referenceBroadcast,
                                 final Broadcast<SAMSequenceDictionary> referenceSequenceDictionaryBroadcast,
                                 final Broadcast<SAMFileHeader> headerBroadcast,
                                 final Broadcast<SVIntervalTree<VariantContext>> cnvCallsBroadcast) {

        Utils.validate(! (evidenceTargetLinks != null && metadata == null),
                "Must supply read metadata when incorporating evidence target links");

        this.sampleId = sampleId;
        this.outputPath = outputPath;
        this.discoverStageArgs = discoverStageArgs;
        this.assemblyRawAlignments = reads;

        this.headerBroadcast = headerBroadcast;
        this.referenceBroadcast = referenceBroadcast;
        this.referenceSequenceDictionaryBroadcast = referenceSequenceDictionaryBroadcast;

        this.cnvCallsBroadcast = cnvCallsBroadcast;
        this.assembledIntervals = assembledIntervals;
        this.evidenceTargetLinks = evidenceTargetLinks;
        this.metadata = metadata;

        this.toolLogger = toolLogger;
    }

    public void updateOutputPath(final String newOutputPath) {
        outputPath = newOutputPath;
    }

    private static Broadcast<SVIntervalTree<VariantContext>> broadcastCNVCalls(final JavaSparkContext ctx,
                                                                               final SAMFileHeader header,
                                                                               final String sampleId,
                                                                               final DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection discoverStageArgs) {
        final SVIntervalTree<VariantContext> cnvCalls;
        if (discoverStageArgs.cnvCallsFile != null) {
            cnvCalls = loadCNVCalls(discoverStageArgs.cnvCallsFile, header, sampleId);
        } else {
            cnvCalls = null;
        }

        final Broadcast<SVIntervalTree<VariantContext>> broadcastCNVCalls;
        if (cnvCalls != null) {
            broadcastCNVCalls = ctx.broadcast(cnvCalls);
        } else {
            broadcastCNVCalls = null;
        }
        return broadcastCNVCalls;
    }

    /**
     * Loads an external cnv call list and returns the results in an SVIntervalTree. NB: the contig indices in the SVIntervalTree
     * are based on the sequence indices in the SAM header, _NOT_ the ReadMetadata (which we might not have access to at this
     * time).
     */
    private static SVIntervalTree<VariantContext> loadCNVCalls(final String cnvCallsFile,
                                                               final SAMFileHeader headerForReads,
                                                               final String sampleId) {
        Utils.validate(cnvCallsFile != null, "Can't load null CNV calls file");
        try ( final FeatureDataSource<VariantContext> dataSource = new FeatureDataSource<>(cnvCallsFile, null, 0, null) ) {
            final VCFHeader cnvCallHeader = (VCFHeader) dataSource.getHeader();
            final ArrayList<String> sampleNamesInOrder = cnvCallHeader.getSampleNamesInOrder();
            Utils.validate(sampleNamesInOrder.size() == 1, "CNV call VCF should be single sample");
            Utils.validate(sampleNamesInOrder.contains(sampleId), ("CNV call VCF does not contain calls for sample " + sampleId));
            Utils.validate(cnvCallHeader.getSequenceDictionary() != null,
                    "CNV calls file does not have a valid sequence dictionary");
            Utils.validate(cnvCallHeader.getSequenceDictionary().isSameDictionary(headerForReads.getSequenceDictionary()),
                    "CNV calls file does not have the same sequence dictionary as the read evidence");
            final SVIntervalTree<VariantContext> cnvCallTree = new SVIntervalTree<>();
            Utils.stream(dataSource.iterator())
                    .map(vc -> new VariantContextBuilder(vc).genotypes(vc.getGenotype(sampleId)).make()) // forces a decode of the genotype for serialization purposes
                    .map(vc -> new Tuple2<>(new SVInterval(headerForReads.getSequenceIndex(vc.getContig()), vc.getStart(), vc.getEnd()),vc))
                    .forEach(pv -> cnvCallTree.put(pv._1(), pv._2()));
            return cnvCallTree;
        }
    }
}
