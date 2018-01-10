package org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.DiscoverVariantsFromContigAlignmentsSAMSpark;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SvCigarUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.io.IOException;
import java.nio.file.Files;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY;

public class AssemblyContigAlignmentsConfigPicker {

    /**
     * Filters an input of SAM file containing alignments of a single-ended long read that
     * aims at providing an "optimal coverage" of the long read, based on an heuristic scoring scheme
     * {@link #computeScoreOfConfiguration(List, Set, int)}.
     *
     * @param longReads             long read alignments
     * @param header                header for the long reads
     * @param scoreDiffTolerance    a tolerance where if two configurations' scores differ less than or equal to this amount, they are considered equally good
     * @param toolLogger            logger for, most likely, debugging uses
     *
     * @return              contigs with alignments filtered and custom formatted as {@link AlignmentInterval}
     */
    public static JavaRDD<AlignedContig> parseIntoOptimallyCoveredContig(final JavaRDD<GATKRead> longReads,
                                                                         final SAMFileHeader header,
                                                                         final String nonCanonicalContigNamesFile,
                                                                         final Double scoreDiffTolerance,
                                                                         final Logger toolLogger) {

        final JavaRDD<AlignedContig> parsedContigAlignments = preprocess(longReads, header, toolLogger);

        return filterAndSplitGappedAI(parsedContigAlignments, nonCanonicalContigNamesFile,
                                      header.getSequenceDictionary(), scoreDiffTolerance);
    }

    private static JavaRDD<AlignedContig> preprocess(final JavaRDD<GATKRead> longReads, final SAMFileHeader header, final Logger toolLogger) {
        longReads.cache();
        toolLogger.info( "Processing " + longReads.count() + " raw alignments from " +
                         longReads.map(GATKRead::getName).distinct().count() + " contigs.");

        final JavaRDD<AlignedContig> parsedContigAlignments =
                new DiscoverVariantsFromContigAlignmentsSAMSpark.SAMFormattedContigAlignmentParser(longReads, header, false)
                        .getAlignedContigs()
                        .filter(AssemblyContigAlignmentsConfigPicker::contigNaiveFilter).cache();
        longReads.unpersist();
        toolLogger.info( "Primitive filtering based purely on MQ left " + parsedContigAlignments.count() + " contigs.");
        return parsedContigAlignments;
    }

    /**
     * Idea is to keep mapped contig that has at least two alignments over MQ 20,
     * or in the case of a single alignment (cannot be simply filtered out the contig because it may contain a large gap),
     * it must be MQ>20.
     */
    private static boolean contigNaiveFilter(final AlignedContig contig) {
        if (contig.alignmentIntervals.size() < 2 ) {
            return (!contig.alignmentIntervals.isEmpty()) && contig.alignmentIntervals.get(0).mapQual > 20;
        } else {
            return contig.alignmentIntervals.stream().mapToInt(ai -> ai.mapQual).filter(mq -> mq > 20).count() > 1;
        }
    }

    /**
     * Manages configuration score-based filtering,
     * then split the gapped alignments by delegating to {@link ContigAlignmentsModifier#splitGappedAlignment(AlignmentInterval, int, int)}.
     * Note that this step is essentially a flatMap operation, meaning that one contig may return 1 or multiple contigs:
     *  *) when 1 contig is yielded, it means the contig has only 1 configuration that scored the best
     *  *) when multiple contigs are yielded, it means the contig has several top-scored configurations
     * How to handle the second scenario can be treated in a separate logic unit.
     */
    @VisibleForTesting
    static JavaRDD<AlignedContig> filterAndSplitGappedAI(final JavaRDD<AlignedContig> parsedContigAlignments,
                                                         final String nonCanonicalContigNamesFile,
                                                         final SAMSequenceDictionary dictionary,
                                                         final Double scoreDiffTolerance) {

        final Set<String> canonicalChromosomes = getCanonicalChromosomes(nonCanonicalContigNamesFile, dictionary);

        return parsedContigAlignments
                .mapToPair(alignedContig -> new Tuple2<>(alignedContig.contigName,
                        new Tuple2<>(alignedContig.contigSequence, pickBestConfigurations(alignedContig, canonicalChromosomes, scoreDiffTolerance))))
                .flatMap(AssemblyContigAlignmentsConfigPicker::reConstructContigFromPickedConfiguration);
    }

    /**
     * Pick the best configurations based on a heuristic scoring scheme implemented in
     * {@link #computeScoreOfConfiguration(List, Set, int)}.
     * @return a 2-D list, where in the case when multiple configurations are equally top-scored, all such configurations are picked up
     */
    @VisibleForTesting
    static List<List<AlignmentInterval>> pickBestConfigurations(final AlignedContig alignedContig,
                                                                final Set<String> canonicalChromosomes,
                                                                final Double scoreDiffTolerance) {

        // group 1: get max aligner score of mappings to canonical chromosomes and speed up in case of too many mappings
        final int maxCanonicalChrAlignerScore = alignedContig.alignmentIntervals.stream()
                .filter(alignmentInterval -> canonicalChromosomes.contains(alignmentInterval.referenceSpan.getContig()))
                .mapToInt(ai -> ai.alnScore).max().orElse(0); // possible that no mapping to canonical chromosomes

        // speed up if number of alignments is too high (>10)
        // if mapped to canonical chromosomes, MQ must be >10; otherwise, must have AS higher than max canonical aligner score
        final List<AlignmentInterval> alignmentIntervals;
        if (alignedContig.alignmentIntervals.size() > 10) {
            alignmentIntervals = alignedContig.alignmentIntervals.stream()
                    .filter(alignmentInterval -> (!canonicalChromosomes.contains(alignmentInterval.referenceSpan.getContig())
                                                                           && alignmentInterval.alnScore > maxCanonicalChrAlignerScore)
                                                 || alignmentInterval.mapQual>10)
                    .collect(Collectors.toList());
        } else {
            alignmentIntervals = alignedContig.alignmentIntervals;
        }

        final int newMaxCanonicalChrAlignerScore = alignmentIntervals.stream()
                .filter(alignmentInterval -> canonicalChromosomes.contains(alignmentInterval.referenceSpan.getContig()))
                .mapToInt(ai -> ai.alnScore).max().orElse(0); // possible that no mapping to canonical chromosomes

        // group 2: generate, and score configurations
        final List<List<AlignmentInterval>> allConfigurations = Sets.powerSet(new HashSet<>(alignmentIntervals))
                .stream().map(ArrayList::new)
                // make sure within each configuration, alignments would be sorted as they would be in a corresponding AlignedContig
                .map(ls -> ls.stream().sorted(AlignedContig.getAlignmentIntervalComparator()).collect(Collectors.toList()))
                .collect(Collectors.toList());

        final List<Double> scores = allConfigurations.stream()
                .map(configuration -> computeScoreOfConfiguration(configuration, canonicalChromosomes, newMaxCanonicalChrAlignerScore))
                .collect(SVUtils.arrayListCollector(allConfigurations.size()));

        // group 3: pick the best-scored configuration(s) (if multiple configurations have equally good scores, return all of them)
        final double maxScore = scores.stream().mapToDouble(Double::doubleValue).max()
                .orElseThrow(() -> new GATKException("Cannot find best-scoring configuration on alignments of contig: " + alignedContig.contigName));

        return IntStream.range(0, allConfigurations.size())
                .filter(i -> {
                    final Double s = scores.get(i);
                    // two configurations with would-be-same-scores can differ by a tolerance due to the sin of comparing
                    // could-be-close floating point values
                    // (see http://www.cygnus-software.com/papers/comparingfloats/Comparing%20floating%20point%20numbers.htm)
                    final Double tol = Math.max(Math.ulp(s), scoreDiffTolerance);
                    return s >= maxScore || maxScore - s <= tol;
                })
                .mapToObj(allConfigurations::get).collect(Collectors.toList());
    }

    /**
     * quick and dirty implementation for computing score of given configuration of alignments,
     * no assumption on the ordering of input alignments
     */
    @VisibleForTesting
    static double computeScoreOfConfiguration(final List<AlignmentInterval> configuration,
                                              final Set<String> canonicalChromosomes,
                                              final int maxCanonicalChrAlignerScore) {

        final double tigExplainQual = computeTigExplainQualOfOneAlignment(configuration, canonicalChromosomes, maxCanonicalChrAlignerScore);

        int redundancy = 0;
        for (int i = 0; i < configuration.size() -1 ; ++i) {
            for (int j = i + 1; j < configuration.size(); ++j) {
                final int overlap = AlignmentInterval.overlapOnContig(configuration.get(i), configuration.get(j));
                redundancy += overlap;
            }
        }

        return tigExplainQual - redundancy;
    }

    private static double computeTigExplainQualOfOneAlignment(final List<AlignmentInterval> configuration,
                                                              final Set<String> canonicalChromosomes,
                                                              final int maxCanonicalChrAlignerScore) {
        double tigExplainedQual = 0;
        for (final AlignmentInterval alignmentInterval : configuration) {
            final int len = alignmentInterval.endInAssembledContig - alignmentInterval.startInAssembledContig + 1;
            final double weight;
            if (canonicalChromosomes.contains(alignmentInterval.referenceSpan.getContig())) {
                weight = alignmentInterval.mapQual/60.0;
            } else {
                weight = Math.max(alignmentInterval.mapQual/60.0,
                                  alignmentInterval.alnScore > maxCanonicalChrAlignerScore ? 1 : 0);
            }
            tigExplainedQual += weight * len;
        }
        return tigExplainedQual;
    }

    /**
     * Depending on the number of best-scored configurations, one contig may produce multiple contigs with the same name
     * and sequence but different selected configurations.
     */
    private static Iterator<AlignedContig> reConstructContigFromPickedConfiguration(
            final Tuple2<String, Tuple2<byte[], List<List<AlignmentInterval>>>> nameSeqAndAlignmentsOfOneRead) {
        if (nameSeqAndAlignmentsOfOneRead._2._2.size() > 1) {
            return nameSeqAndAlignmentsOfOneRead._2._2.stream()
                    .map(configuration -> new AlignedContig(nameSeqAndAlignmentsOfOneRead._1, nameSeqAndAlignmentsOfOneRead._2._1,
                            splitGaps(configuration), true))
                    .sorted(sortConfigurations())
                    .collect(Collectors.toList()).iterator();
        } else {
            return Collections.singletonList(new AlignedContig(nameSeqAndAlignmentsOfOneRead._1, nameSeqAndAlignmentsOfOneRead._2._1,
                    splitGaps(nameSeqAndAlignmentsOfOneRead._2._2.get(0)), false))
                    .iterator();
        }
    }

    /**
     * when two configurations are the same, prefer the one with less alignments,
     * or less summed mismatches if still tie.
     */
    private static Comparator<AlignedContig> sortConfigurations() {
        final Comparator<AlignedContig> numFirst
                = (AlignedContig x, AlignedContig y) -> Integer.compare(x.alignmentIntervals.size(), y.alignmentIntervals.size());
        final Comparator<AlignedContig> mismatchSecond
                = (AlignedContig x, AlignedContig y) -> Integer.compare(x.alignmentIntervals.stream().mapToInt(ai -> ai.mismatches).sum(),
                                                                        y.alignmentIntervals.stream().mapToInt(ai -> ai.mismatches).sum());
        return numFirst.thenComparing(mismatchSecond);
    }

    private static List<AlignmentInterval> splitGaps(final List<AlignmentInterval> configuration) {
        return configuration.stream()
                .map(ai -> ContigAlignmentsModifier.splitGappedAlignment(ai, GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY,
                        SvCigarUtils.getUnclippedReadLength(ai.cigarAlong5to3DirectionOfContig)))
                .flatMap(Utils::stream).collect(Collectors.toList());
    }

    /**
     * Primary reference contigs are defined as chromosomes 1-22, X, Y, M, and defined for both GRCh38 and hg19.
     */
    static Set<String> getCanonicalChromosomes(final String nonCanonicalContigNamesFile, final SAMSequenceDictionary dictionary) {
        if (nonCanonicalContigNamesFile!= null) {

            try (final Stream<String> nonCanonical = Files.lines(IOUtils.getPath((nonCanonicalContigNamesFile)))) {
                return new HashSet<>( Sets.difference(dictionary.getSequences().stream().map(SAMSequenceRecord::getSequenceName).collect(Collectors.toSet()),
                        nonCanonical.collect(Collectors.toSet())) );
            } catch ( final IOException ioe ) {
                throw new UserException("Can't read nonCanonicalContigNamesFile file "+nonCanonicalContigNamesFile, ioe);
            }
        } else {
            final List<String> first22ChromosomesNum = IntStream.range(0, 23).mapToObj(String::valueOf).collect(Collectors.toList());
            final Set<String> canonicalChromosomeNames = first22ChromosomesNum.stream().map(name -> "chr" + name).collect(Collectors.toSet());
            canonicalChromosomeNames.addAll(first22ChromosomesNum);
            canonicalChromosomeNames.addAll(Arrays.asList("chrX", "chrY", "chrM", "X", "Y", "MT"));
            return new HashSet<>( canonicalChromosomeNames );
        }
    }
}
