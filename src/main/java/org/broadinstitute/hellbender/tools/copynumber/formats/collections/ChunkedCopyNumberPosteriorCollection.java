package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import org.broadinstitute.hellbender.tools.copynumber.gcnv.CopyNumberPosteriorRecord;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.IntegerCopyNumberState;
import org.broadinstitute.hellbender.tools.copynumber.gcnv.IntegerCopyNumberStateCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.BiConsumer;
import java.util.function.Function;

/**
 * Collection of copy number posteriors for an individual chunk containing a subset of intervals considered in the analysis
 */
public class ChunkedCopyNumberPosteriorCollection extends AbstractSampleLocatableCollection<CopyNumberPosteriorRecord> {

    public ChunkedCopyNumberPosteriorCollection(final File inputFile,
                                                final List<SimpleInterval> chunkIntervals,
                                                final IntegerCopyNumberStateCollection integerCopyNumberStateCollection) {
        super(inputFile,
                Utils.nonNull(integerCopyNumberStateCollection).getTableColumnCollection(),
                getPosteriorRecordFromDataLineDecoder(chunkIntervals, integerCopyNumberStateCollection),
                getPosteriorRecordToDataLineEncoder(chunkIntervals, integerCopyNumberStateCollection));
    }

    private static Function<DataLine, CopyNumberPosteriorRecord> getPosteriorRecordFromDataLineDecoder(final List<SimpleInterval> chunkIntervals,
                                                                                                       final IntegerCopyNumberStateCollection integerCopyNumberStateCollection) {
        return dataLine -> {
            final SimpleInterval interval = chunkIntervals.get((int) dataLine.getLineNumber());
            final Map<IntegerCopyNumberState, Double> copyNumberStateDoubleMap = new HashMap<>();
            for(int i = 0; i < integerCopyNumberStateCollection.size(); i++) {
                copyNumberStateDoubleMap.putIfAbsent(integerCopyNumberStateCollection.get(i), dataLine.getDouble(i));
            }
            final CopyNumberPosteriorRecord record = new CopyNumberPosteriorRecord(interval, copyNumberStateDoubleMap);
            return record;
        };
    }

    private static BiConsumer<CopyNumberPosteriorRecord, DataLine> getPosteriorRecordToDataLineEncoder(final List<SimpleInterval> chunkIntervals,
                                                                                                       final IntegerCopyNumberStateCollection integerCopyNumberStateCollection) {
        //TODO
        return (copyNumberPosteriorRecord, dataLine) -> {
            return;
        };
    }

}
