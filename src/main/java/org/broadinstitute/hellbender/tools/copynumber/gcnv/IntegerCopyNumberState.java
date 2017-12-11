package org.broadinstitute.hellbender.tools.copynumber.gcnv;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.param.ParamUtils;


/**
 * This class represents integer copy number states. It is intended to be used as states of a
 * hidden Markov model.
 */
public final class IntegerCopyNumberState  {

    /**
     * Integer value of the represented copy number state
     */
    private final int copyNumber;

    /**
     * A string representation of this copy number state (used for creating human-readable copy number call lists)
     */
    private final String callString;

    /**
     * An allele representation of this copy number state (used for VCF creation)
     */
    private final Allele allele;

    public IntegerCopyNumberState(final int copyNumber) {
        this.copyNumber = ParamUtils.isPositiveOrZero(copyNumber, "The integer copy number state" +
                " must be non-negative");
        this.callString = toCallString(copyNumber);
        allele = Allele.create(toAlleleString(copyNumber));
    }

    public Allele toAllele() {
        return allele;
    }

    public String getCallString() {
        return callString;
    }

    public double getScalar() {
        return copyNumber;
    }

    public int getCopyNumber() { return copyNumber; }

    private static String toCallString(final int copyNumberState) {
        return String.valueOf(copyNumberState);
    }

    private static String toAlleleString(final int copyNumberState) {
        return "<" + String.valueOf(copyNumberState) + ">";
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        return copyNumber == ((IntegerCopyNumberState) o).copyNumber;
    }

    @Override
    public int hashCode() {
        return copyNumber;
    }
}