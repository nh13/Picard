package picard.metrics;


import picard.analysis.replicates.MergeableMetricBase;

public class MultilevelMetrics extends MergeableMetricBase {
     /** The sample to which these metrics apply.  If null, it means they apply
     * to all reads in the file. */
     @MergeByAssertEquals
    public String SAMPLE;

    /** The library to which these metrics apply.  If null, it means that the
     * metrics were accumulated at the sample level. */
    @MergeByAssertEquals
    public String LIBRARY = null;

    /** The read group to which these metrics apply.  If null, it means that
     * the metrics were accumulated at the library or sample level.*/
    @MergeByAssertEquals
    public String READ_GROUP = null;
}
