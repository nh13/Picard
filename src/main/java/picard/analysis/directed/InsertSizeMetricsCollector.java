package picard.analysis.directed;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamPairUtil;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Histogram;
import picard.analysis.InsertSizeMetrics;
import picard.analysis.MetricAccumulationLevel;
import picard.metrics.MultiLevelCollector;
import picard.metrics.PerUnitMetricCollector;

import java.util.Collections;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Collects InsertSizeMetrics on the specified accumulationLevels using
 */
public class InsertSizeMetricsCollector extends MultiLevelCollector<InsertSizeMetrics, Integer, InsertSizeCollectorArgs> {
    // When generating the Histogram, discard any data categories (out of FR, TANDEM, RF) that have fewer than this
    // percentage of overall reads. (Range: 0 to 1)
    private final double minimumPct;

    // Generate mean, sd and plots by trimming the data down to MEDIAN + DEVIATIONS*MEDIAN_ABSOLUTE_DEVIATION.
    // This is done because insert size data typically includes enough anomalous values from chimeras and other
    // artifacts to make the mean and sd grossly misleading regarding the real distribution.
    private final double deviations;

    //Explicitly sets the Histogram width, overriding automatic truncation of Histogram tail.
    //Also, when calculating mean and stdev, only bins <= Histogram_WIDTH will be included.
    private final Integer histogramWidth;

    // If set to true, then duplicates will also be included in the histogram
    private final boolean includeDuplicates;

    public InsertSizeMetricsCollector(final Set<MetricAccumulationLevel> accumulationLevels, final List<SAMReadGroupRecord> samRgRecords,
                                      final double minimumPct, final Integer histogramWidth, final double deviations,
                                      final boolean includeDuplicates) {
        this.minimumPct = minimumPct;
        this.histogramWidth = histogramWidth;
        this.deviations = deviations;
        this.includeDuplicates = includeDuplicates;
        setup(accumulationLevels, samRgRecords);
    }

    // We will pass insertSize and PairOrientation with the DefaultPerRecordCollectorArgs passed to the record collectors
    // This method is called once Per samRecord
    @Override
    protected InsertSizeCollectorArgs makeArg(SAMRecord samRecord, ReferenceSequence refSeq) {
        final int insertSize = Math.abs(samRecord.getInferredInsertSize());
        final SamPairUtil.PairOrientation orientation = SamPairUtil.getPairOrientation(samRecord);

        return new InsertSizeCollectorArgs(insertSize, orientation);
    }

    /** Make an InsertSizeCollector with the given arguments */
    @Override
    protected PerUnitMetricCollector<InsertSizeMetrics, Integer, InsertSizeCollectorArgs> makeChildCollector(final String sample, final String library, final String readGroup) {
        return new PerUnitInsertSizeMetricsCollector(sample, library, readGroup);
    }

    @Override
    public void acceptRecord(final SAMRecord record, final ReferenceSequence refSeq) {
        if (!record.getReadPairedFlag() ||
                record.getReadUnmappedFlag() ||
                record.getMateUnmappedFlag() ||
                record.getFirstOfPairFlag() ||
                record.isSecondaryOrSupplementary() ||
                (record.getDuplicateReadFlag() && !this.includeDuplicates) ||
                record.getInferredInsertSize() == 0) {
            return;
        }

        super.acceptRecord(record, refSeq);
    }

    /** A Collector for individual InsertSizeMetrics for a given SAMPLE or SAMPLE/LIBRARY or SAMPLE/LIBRARY/READ_GROUP (depending on aggregation levels) */
    public class PerUnitInsertSizeMetricsCollector implements PerUnitMetricCollector<InsertSizeMetrics, Integer, InsertSizeCollectorArgs> {
        final EnumMap<SamPairUtil.PairOrientation, Histogram<Integer>> histograms = new EnumMap<SamPairUtil.PairOrientation, Histogram<Integer>>(SamPairUtil.PairOrientation.class);
        final String sample;
        final String library;
        final String readGroup;
        private double totalInserts = 0;

        public PerUnitInsertSizeMetricsCollector(final String sample, final String library, final String readGroup) {
            this.sample = sample;
            this.library = library;
            this.readGroup = readGroup;
            String prefix = null;
            if (this.readGroup != null) {
                prefix = this.readGroup + ".";
            }
            else if (this.library != null) {
                prefix = this.library + ".";
            }
            else if (this.sample != null) {
                prefix = this.sample + ".";
            }
            else {
                prefix = "All_Reads.";
            }
            histograms.put(SamPairUtil.PairOrientation.FR,     new Histogram<Integer>("insert_size", prefix + "fr_count"));
            histograms.put(SamPairUtil.PairOrientation.TANDEM, new Histogram<Integer>("insert_size", prefix + "tandem_count"));
            histograms.put(SamPairUtil.PairOrientation.RF,     new Histogram<Integer>("insert_size", prefix + "rf_count"));
        }

        public void acceptRecord(final InsertSizeCollectorArgs args) {
            histograms.get(args.getPairOrientation()).increment(args.getInsertSize());
        }

        public void finish() { }

        public double getTotalInserts() {
            return totalInserts;
        }

        public void addMetricsToFile(final MetricsFile<InsertSizeMetrics,Integer> file) {
            // get the number of inserts, and the maximum and minimum keys across, across all orientations
            for (final Histogram<Integer> h : this.histograms.values()) {
                totalInserts += h.getCount();
            }
            if (0 == totalInserts) return; // nothing to store

            for(final Map.Entry<SamPairUtil.PairOrientation, Histogram<Integer>> entry : histograms.entrySet()) {
                final SamPairUtil.PairOrientation pairOrientation = entry.getKey();
                final Histogram<Integer> histogram = entry.getValue();
                final double total = histogram.getCount();

                // Only include a category if it has a sufficient percentage of the data in it
                if( total >= totalInserts * minimumPct ) {
                    final InsertSizeMetrics metrics = new InsertSizeMetrics();
                    metrics.SAMPLE             = this.sample;
                    metrics.LIBRARY            = this.library;
                    metrics.READ_GROUP         = this.readGroup;
                    metrics.PAIR_ORIENTATION   = pairOrientation;

                    if (!histogram.isEmpty()) {
                        metrics.READ_PAIRS = (long) total;
                        metrics.calculateDerivedFields(histogram, getWidthToTrimTo(metrics));
                    }

                    file.addHistogram(histogram);
                    file.addMetric(metrics);
                }
            }
        }

        /**
         * @return {@link #histogramWidth} if it was specified in the constructor or a calculated width based on the stdev of the input metric and {@link #deviations}
         */
        private int getWidthToTrimTo(InsertSizeMetrics metrics) {
            if (histogramWidth == null) {
                return metrics.getWidthToTrimTo(deviations);
            } else {
                return histogramWidth;
            }
        }
    }
}

// Arguments that need to be calculated once per SAMRecord that are then passed to each PerUnitMetricCollector
// for the given record
class InsertSizeCollectorArgs {
    private final int insertSize;
    private final SamPairUtil.PairOrientation po;


    public int getInsertSize() {
        return insertSize;
    }

    public SamPairUtil.PairOrientation getPairOrientation() {
        return po;
    }

    public InsertSizeCollectorArgs(final int insertSize, final SamPairUtil.PairOrientation po) {
        this.insertSize = insertSize;
        this.po = po;
    }
}
