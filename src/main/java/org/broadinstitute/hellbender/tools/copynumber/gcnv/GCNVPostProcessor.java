package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.hellbender.engine.ProgressMeter;
import org.broadinstitute.hellbender.utils.Utils;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Optional;

/**
 * Helper class for single sample gCNV postprocessing
 */
public class GCNVPostProcessor {

    /**
     * VCF header keys
     */
    public static final String CN_MAP = "CNMAP";
    public static final String CN_MEAN = "CNMEAN";
    public static final String CN_STD = "CNSTD";
    public static final String CNQ = "CNQ";

    private final List<Allele> alleles;
    private final IntegerCopyNumberStateCollection integerCopyNumberStateCollection;
    private final String sampleName;
    private final VariantContextWriter outputWriter;

    GCNVPostProcessor(final VariantContextWriter outputWriter,
                      final IntegerCopyNumberStateCollection integerCopyNumberStateCollection,
                      final String sampleName) {
        this.outputWriter = Utils.nonNull(outputWriter);
        this.integerCopyNumberStateCollection = Utils.nonNull(integerCopyNumberStateCollection);
        this.alleles = integerCopyNumberStateCollection.getAlleles();
        this.sampleName = sampleName;
    }

    /**
     *
     */
    public void composeVariantContextAndWrite(@Nullable final String commandLine) {
        //TODO pass a list of intervals and add progress meter
        //ProgressMeter progressMeter = new ProgressMeter(1.0);
        //progressMeter.start();
        outputWriter.writeHeader(composeHeader(commandLine));
    }

    public void writeChunkedVariantContext(final List<CopyNumberPosteriorRecord> copyNumberPosteriorRecordsChunk,
                                           final String variantPrefix) {
        for (CopyNumberPosteriorRecord posteriorRecord: copyNumberPosteriorRecordsChunk) {
            final VariantContext variantContext = composeVariantContext(posteriorRecord, variantPrefix);
            outputWriter.add(variantContext);
        }
    }

    /**
     * Composes the VCF header
     */
    private VCFHeader composeHeader(@Nullable final String commandLine) {
        final VCFHeader result = new VCFHeader(Collections.emptySet(), Arrays.asList(sampleName));

        /* add VCF version */
        result.addMetaDataLine(new VCFHeaderLine(VCFHeaderVersion.VCF4_2.getFormatString(),
                VCFHeaderVersion.VCF4_2.getVersionString()));

        /* add command line */
        if (commandLine != null) {
            result.addMetaDataLine(new VCFHeaderLine("command", commandLine));
        }

        /* header lines related to formatting */
        result.addMetaDataLine(new VCFFormatHeaderLine(CN_MAP, 1,
                VCFHeaderLineType.Integer, "Copy number maximum a posteriori"));
        result.addMetaDataLine(new VCFFormatHeaderLine(CN_MEAN, 1,
                VCFHeaderLineType.Float, "Copy number posterior mean"));

        return result;
    }

    /**
     *
     * @param copyNumberPosteriorRecord
     * @param variantPrefix
     * @return
     */
    protected VariantContext composeVariantContext(final CopyNumberPosteriorRecord copyNumberPosteriorRecord,
                                                 final String variantPrefix) {
        final VariantContextBuilder variantContextBuilder = new VariantContextBuilder();
        variantContextBuilder.alleles(alleles);
        variantContextBuilder.chr(copyNumberPosteriorRecord.getContig());
        variantContextBuilder.start(copyNumberPosteriorRecord.getStart());
        variantContextBuilder.stop(copyNumberPosteriorRecord.getEnd());
        variantContextBuilder.id(String.format(variantPrefix + "_%s_%d_%d",
                copyNumberPosteriorRecord.getContig(),
                copyNumberPosteriorRecord.getStart(),
                copyNumberPosteriorRecord.getEnd()));
        final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sampleName);
        final int copyNumberMAP = calculateMAPCopyNumber(copyNumberPosteriorRecord);
        genotypeBuilder.attribute(CN_MAP, copyNumberMAP);
        final Genotype genotype = genotypeBuilder.make();
        //TODO add variant context fields to genotype

        //Add allele information to the variant context
        variantContextBuilder.attribute(VCFConstants.END_KEY, copyNumberPosteriorRecord.getEnd());

        variantContextBuilder.genotypes(genotype);
        return variantContextBuilder.make();
    }

    private int calculateMAPCopyNumber(CopyNumberPosteriorRecord copyNumberPosteriorRecord) {
        Optional<IntegerCopyNumberState> copyNumberStateMAP = integerCopyNumberStateCollection.getCopyNumberStates()
                .stream().max((a1, a2) -> Double.compare(copyNumberPosteriorRecord.getCopyNumberPosterior(a1),
                copyNumberPosteriorRecord.getCopyNumberPosterior(a2)));
        return copyNumberStateMAP.get().getCopyNumber();
    }


}
