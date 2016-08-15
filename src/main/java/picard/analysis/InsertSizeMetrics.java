/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package picard.analysis;

import htsjdk.samtools.SamPairUtil.PairOrientation;
import htsjdk.samtools.util.Histogram;
import picard.PicardException;
import picard.metrics.MultilevelMetrics;

/**
 * Metrics about the insert size distribution of a paired-end library, created by the
 * CollectInsertSizeMetrics program and usually written to a file with the extension
 * ".insert_size_metrics".  In addition the insert size distribution is plotted to
 * a file with the extension ".insert_size_Histogram.pdf".
 *
 * @author Doug Voet (dvoet at broadinstitute dot org)
 */
public class InsertSizeMetrics extends MultilevelMetrics {

    /** The MEDIAN insert size of all paired end reads where both ends mapped to the same chromosome. */
    @NoMergingIsDerived
    public double MEDIAN_INSERT_SIZE;

    /**
     * The median absolute deviation of the distribution.  If the distribution is essentially normal then
     * the standard deviation can be estimated as ~1.4826 * MAD.
     */
    @NoMergingIsDerived
    public double MEDIAN_ABSOLUTE_DEVIATION;

    /** The minimum measured insert size.  This is usually 1 and not very useful as it is likely artifactual. */
    @NoMergingIsDerived
    public int MIN_INSERT_SIZE;
    /**
     * The maximum measure insert size by alignment. This is usually very high representing either an artifact
     * or possibly the presence of a structural re-arrangement.
     */
    @NoMergingIsDerived
    public int MAX_INSERT_SIZE;
    /**
     * The mean insert size of the "core" of the distribution. Artefactual outliers in the distribution often
     * cause calculation of nonsensical mean and stdev values.  To avoid this the distribution is first trimmed
     * to a "core" distribution of +/- N median absolute deviations around the median insert size. By default
     * N=10, but this is configurable.
     */
    @NoMergingIsDerived
    public double MEAN_INSERT_SIZE;
    /** Standard deviation of insert sizes over the "core" of the distribution. */
    @NoMergingIsDerived
    public double STANDARD_DEVIATION;
    /** The total number of read pairs that were examined in the entire distribution. */
    @MergeByAdding
    public long READ_PAIRS;
    /** The pair orientation of the reads in this data category. */
    @MergeByAssertEquals
    public PairOrientation PAIR_ORIENTATION;

    /** The "width" of the bins, centered around the median, that encompass 10% of all read pairs. */
    @NoMergingIsDerived
    public int WIDTH_OF_10_PERCENT;
    /** The "width" of the bins, centered around the median, that encompass 20% of all read pairs. */
    @NoMergingIsDerived
    public int WIDTH_OF_20_PERCENT;
    /** The "width" of the bins, centered around the median, that encompass 30% of all read pairs. */
    @NoMergingIsDerived
    public int WIDTH_OF_30_PERCENT;
    /** The "width" of the bins, centered around the median, that encompass 40% of all read pairs. */
    @NoMergingIsDerived
    public int WIDTH_OF_40_PERCENT;
    /** The "width" of the bins, centered around the median, that encompass 50% of all read pairs. */
    @NoMergingIsDerived
    public int WIDTH_OF_50_PERCENT;
    /** The "width" of the bins, centered around the median, that encompass 60% of all read pairs. */
    @NoMergingIsDerived
    public int WIDTH_OF_60_PERCENT;
    /**
     * The "width" of the bins, centered around the median, that encompass 70% of all read pairs.
     * This metric divided by 2 should approximate the standard deviation when the insert size
     * distribution is a normal distribution.
     */
    @NoMergingIsDerived
    public int WIDTH_OF_70_PERCENT;
    /** The "width" of the bins, centered around the median, that encompass 80% of all read pairs. */
    @NoMergingIsDerived
    public int WIDTH_OF_80_PERCENT;
    /** The "width" of the bins, centered around the median, that encompass 90% of all read pairs. */
    @NoMergingIsDerived
    public int WIDTH_OF_90_PERCENT;
    /** The "width" of the bins, centered around the median, that encompass 100% of all read pairs. */
    @NoMergingIsDerived
    public int WIDTH_OF_99_PERCENT;


    /** Gets the width of the histogram for trimming based on the stdev of the input metric. */
    public int getWidthToTrimTo(final double deviations) {
        return (int) (this.MEDIAN_INSERT_SIZE + (deviations * this.MEDIAN_ABSOLUTE_DEVIATION));
    }

    /**
     * Not implemented, use {@link #calculateDerivedFields(Histogram, int)} instead.
     */
    @Override
    public void calculateDerivedFields() {
        throw new PicardException("Use calculateDerivedFields(");
    }


    /**
     * Calculate the derived fields for this metric.
     *
     * @param histogram the histogram for the insert sizes. This histogram should not be trimmed.
     * @param histogramWidth the width of the histogram centered around the median to trim to when calculating the mean
     *                       and standard deviation. This is done because insert size data typically includes enough
     *                       anomalous values from chimeras and other artifacts to make the mean and sd grossly
     *                       misleading regarding the real distribution.
     */
    public void calculateDerivedFields(final Histogram<Integer> histogram, final int histogramWidth) {

        if (!histogram.isEmpty()) {
            MAX_INSERT_SIZE = (int) histogram.getMax();
            MIN_INSERT_SIZE = (int) histogram.getMin();
            MEDIAN_INSERT_SIZE = histogram.getMedian();
            MEDIAN_ABSOLUTE_DEVIATION = histogram.getMedianAbsoluteDeviation();

            final double median = histogram.getMedian();
            double covered = 0;
            double low = median;
            double high = median;

            while (low >= histogram.getMin() || high <= histogram.getMax()) {
                final Histogram.Bin<Integer> lowBin = histogram.get((int) low);
                if (lowBin != null) covered += lowBin.getValue();

                if (low != high) {
                    final Histogram.Bin<Integer> highBin = histogram.get((int) high);
                    if (highBin != null) covered += highBin.getValue();
                }

                final double percentCovered = covered / READ_PAIRS;
                final int distance = (int) (high - low) + 1;
                if (percentCovered >= 0.1 && WIDTH_OF_10_PERCENT == 0) WIDTH_OF_10_PERCENT = distance;
                if (percentCovered >= 0.2 && WIDTH_OF_20_PERCENT == 0) WIDTH_OF_20_PERCENT = distance;
                if (percentCovered >= 0.3 && WIDTH_OF_30_PERCENT == 0) WIDTH_OF_30_PERCENT = distance;
                if (percentCovered >= 0.4 && WIDTH_OF_40_PERCENT == 0) WIDTH_OF_40_PERCENT = distance;
                if (percentCovered >= 0.5 && WIDTH_OF_50_PERCENT == 0) WIDTH_OF_50_PERCENT = distance;
                if (percentCovered >= 0.6 && WIDTH_OF_60_PERCENT == 0) WIDTH_OF_60_PERCENT = distance;
                if (percentCovered >= 0.7 && WIDTH_OF_70_PERCENT == 0) WIDTH_OF_70_PERCENT = distance;
                if (percentCovered >= 0.8 && WIDTH_OF_80_PERCENT == 0) WIDTH_OF_80_PERCENT = distance;
                if (percentCovered >= 0.9 && WIDTH_OF_90_PERCENT == 0) WIDTH_OF_90_PERCENT = distance;
                if (percentCovered >= 0.99 && WIDTH_OF_99_PERCENT == 0) WIDTH_OF_99_PERCENT = distance;

                --low;
                ++high;
            }
        }

        // Trim the Histogram down to get rid of outliers that would make the mean and stddev useless.
        histogram.trimByWidth(histogramWidth);

        if (!histogram.isEmpty()) {
            MEAN_INSERT_SIZE = histogram.getMean();
            STANDARD_DEVIATION = histogram.getStandardDeviation();
        }
    }
}
