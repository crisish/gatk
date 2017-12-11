package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.ChunkedCopyNumberPosteriorCollection;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.text.XReadLines;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

/**
 * This tool takes in a list of directories to chunked gCNV calls output containing posteriors for different copy number states,
 * and outputs a VCF containing records
 */
@CommandLineProgramProperties(
        summary = "Collects ref/alt counts at sites.",
        oneLineSummary = "Collects ref/alt counts at sites.",
        programGroup = CopyNumberProgramGroup.class
)
public class GenerateVCFFromPosteriors extends GATKTool {
    private static final Logger logger = LogManager.getLogger(GenerateVCFFromPosteriors.class);

    public static final String POSTERIOR_CALL_DIRECTORY_FULL_NAME = "copy-number-posterior";
    public static final String POSTERIOR_CALL_DIRECTORY_SHORT_NAME = "CNP";

    public static final String SAMPLE_NAME_FULL_NAME = "sample-name";
    public static final String SAMPLE_NAME_SHORT_NAME = "SN";

    @Argument(
            doc = "List of paths to chunks' directories",
            fullName = POSTERIOR_CALL_DIRECTORY_FULL_NAME,
            shortName = POSTERIOR_CALL_DIRECTORY_SHORT_NAME
    )
    private List<String> chunkDirectoryList;

    @Argument(
            doc = "Index of the sample in the gCNV pipeline",
            fullName = SAMPLE_NAME_FULL_NAME,
            shortName = SAMPLE_NAME_SHORT_NAME
    )
    private String sampleName;

    @Argument(
            doc = "Output VCF file",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputFile;

    /**
     * Genomic location parser
     */
    private GenomeLocParser genomeLocParser;

    final String variantPrefix = "CNV";

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void traverse() {
        logger.info("Initializing and validating intervals...");
        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        genomeLocParser = new GenomeLocParser(sequenceDictionary);
        final VariantContextWriter outputWriter = GATKVariantContextUtils.createVCFWriter(outputFile, null, false);

        //final List<SimpleInterval> intervals = intervalArgumentCollection.getIntervals(sequenceDictionary);

        if (chunkDirectoryList.size() <= 0) {
            throw new UserException.BadInput("You must specify at least one chunk directory");
        }

        Pair<ChunkedCopyNumberPosteriorCollection, IntegerCopyNumberStateCollection> firstCopyNumberPosteriorFileInfo =
                readChunkedPosteriorFileFromDirectory(chunkDirectoryList.get(0));
        final String sampleName = firstCopyNumberPosteriorFileInfo.getLeft().getMetadata().getSampleName(); //TODO check that metadata has the right type
        final IntegerCopyNumberStateCollection copyNumberStateCollection = firstCopyNumberPosteriorFileInfo.getRight();

        final GCNVPostProcessor gcnvPostProcessor = new GCNVPostProcessor(outputWriter, copyNumberStateCollection, sampleName);
        gcnvPostProcessor.composeVariantContextAndWrite(this.getCommandLine());
        int currentChunk = 0;

        for(String chunkRootDirectory: chunkDirectoryList) {
            logger.info(String.format("Analyzing copy number posterior chunk number %d", currentChunk));

            Pair<ChunkedCopyNumberPosteriorCollection, IntegerCopyNumberStateCollection> chunkedCopyNumberPosteriorFileInfo =
            readChunkedPosteriorFileFromDirectory(chunkRootDirectory);
            //check that the sample name and copy number state collection is equivalent across all chunks
            if (!sampleName.equals(chunkedCopyNumberPosteriorFileInfo.getLeft().getSampleName())) {
                throw new UserException.BadInput("Sample name differs across chunked posterior outputs");
            }
            if (!copyNumberStateCollection.equals(chunkedCopyNumberPosteriorFileInfo.getRight())) {
                throw new UserException.BadInput("Copy number collection differs across chunked posterior outputs");
            }

            final List<CopyNumberPosteriorRecord> copyNumberPosteriorRecordsList =
                    chunkedCopyNumberPosteriorFileInfo.getLeft().getRecords();

            gcnvPostProcessor.writeChunkedVariantContext(copyNumberPosteriorRecordsList, variantPrefix);
            currentChunk++;

        }
        return;
    }

    /**
     * Process a single posterior chunked file
     */
    private Pair<ChunkedCopyNumberPosteriorCollection, IntegerCopyNumberStateCollection> readChunkedPosteriorFileFromDirectory(final String chunkRootDirectory) {
        //get a list of intervals associated with chunk currently being processed
        final List<GenomeLoc> chunkGenomeLocList = IntervalUtils.intervalFileToList(genomeLocParser,
                getIntervalFileFromChunkDirectory(chunkRootDirectory).getAbsolutePath());
        final List<SimpleInterval> chunkIntervalList = IntervalUtils.convertGenomeLocsToSimpleIntervals(chunkGenomeLocList);

        final File chunkPosteriorFile = getPosteriorFileFromChunkDirectory(chunkRootDirectory, sampleName);

        final List<String> copyNumberStateColumns = getPosteriorFileColumns(chunkPosteriorFile);
        final IntegerCopyNumberStateCollection chunkCopyNumberStateCollection = new IntegerCopyNumberStateCollection(copyNumberStateColumns);

        final ChunkedCopyNumberPosteriorCollection chunkedCopyNumberPosteriorCollection =
                new ChunkedCopyNumberPosteriorCollection(chunkPosteriorFile, chunkIntervalList, chunkCopyNumberStateCollection);

        return new ImmutablePair<>(chunkedCopyNumberPosteriorCollection, chunkCopyNumberStateCollection);
    }

    /**
     * Get list of column names of the copy number posterior file
     */
    private static List<String> getPosteriorFileColumns(final File copyNumberPosteriorFile) {
        List<String> columns = null;
        try (final XReadLines reader = new XReadLines(copyNumberPosteriorFile)) {
            while (reader.hasNext()) {
                String nextLine = reader.next();
                if (!nextLine.startsWith(TableUtils.COMMENT_PREFIX)) {
                    columns = Arrays.asList(nextLine.split(TableUtils.COLUMN_SEPARATOR_STRING));
                    break;
                }
            }
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(copyNumberPosteriorFile);
        }
        if (columns == null) {
            throw new UserException.BadInput("Copy number posterior file does not have a header");
        }
        return columns;
    }

    /**
     * Get the posterior file
     * TODO extract this to the GermlineCNVNamingConstants class
     */
    private static File getPosteriorFileFromChunkDirectory(final String chunkDirectoryPath, final String sampleIndex) {
        final String posteriorCallsDirectory = chunkDirectoryPath + File.separator + sampleIndex;
        final File posteriorFile = new File(posteriorCallsDirectory, GermlineCNVNamingConstants.COPY_NUMBER_POSTERIOR_FILE_NAME);
        return posteriorFile;
    }

    /**
     * Get the intervals file
     * TODO extract this to the GermlineCNVNamingConstants class
     */
    private static File getIntervalFileFromChunkDirectory(final String chunkDirectoryPath) {
        final File intervalsFile = new File(chunkDirectoryPath, GermlineCNVNamingConstants.INTERVAL_LIST_FILE_NAME);
        return intervalsFile;
    }
}
