package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Map;

/**
 * A single copy number posterior record for a specific interval
 */
public class CopyNumberPosteriorRecord implements Locatable {
    private final SimpleInterval interval;
    private final Map<IntegerCopyNumberState, Double> copyNumberStatePosteriors;

    public CopyNumberPosteriorRecord(final SimpleInterval interval, final Map<IntegerCopyNumberState, Double> copyNumberStatePosteriors) {
        this.interval = Utils.nonNull(interval);
        this.copyNumberStatePosteriors = Utils.nonNull(copyNumberStatePosteriors);
    }

    /**
     * Get the posterior probability for a given copy number state
     */
    public Double getCopyNumberPosterior(IntegerCopyNumberState integerCopyNumberState) {
        return copyNumberStatePosteriors.get(integerCopyNumberState);
    }

    @Override
    public String getContig() {
        return interval.getContig();
    }

    @Override
    public int getStart() {
        return interval.getStart();
    }

    @Override
    public int getEnd() {
        return interval.getEnd();
    }
}
