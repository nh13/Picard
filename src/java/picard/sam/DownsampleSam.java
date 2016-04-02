/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
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
package picard.sam;

import htsjdk.samtools.AbstractBAMFileIndex;
import htsjdk.samtools.BAMIndexMetaData;
import htsjdk.samtools.DownsamplingIterator;
import htsjdk.samtools.DownsamplingIteratorFactory;
import htsjdk.samtools.DownsamplingIteratorFactory.Strategy;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.SamOrBam;

import java.io.File;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.stream.StreamSupport;

/**
 * Class to randomly downsample a BAM file while respecting that we should either retain or discard
 * all of the reads for a template - i.e. all reads with the same name, whether first or second of
 * pair, secondary or supplementary, all travel together.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        usage = DownsampleSam.USAGE_SUMMARY + DownsampleSam.USAGE_DETAILS,
        usageShort = DownsampleSam.USAGE_SUMMARY,
        programGroup = SamOrBam.class
)
public class DownsampleSam extends CommandLineProgram {
    static final String USAGE_SUMMARY = "Downsample a SAM or BAM file.  ";
    static final String USAGE_DETAILS = "This tool applies a random downsampling algorithm to a SAM or BAM file to retain " +
            "only a random subset of the reads. Reads in a mate-pair are either both kept or both discarded. Reads marked as not primary " +
            "alignments are all discarded. Each read is given a probability P of being retained so that runs performed with the exact " +
            "same input in the same order and with the same value for RANDOM_SEED will produce the same results." +
            "All reads for a template are kept or discarded as a unit, with the goal of retaining reads" +
            "from PROBABILITY * input templates. While this will usually result in approximately " +
            "PROBABILITY * input reads being retained also, for very small PROBABILITIES this may not " +
            "be the case.\n" +
            "A number of different downsampling strategies are supported using the STRATEGY option:\n\n" +
            "ConstantMemory: " + DownsamplingIteratorFactory.CONSTANT_MEMORY_DESCRPTION + "\n\n" +
            "HighAccuracy: " + DownsamplingIteratorFactory.HIGH_ACCURACY_DESCRIPTION + "\n\n" +
            "Chained: " + DownsamplingIteratorFactory.CHAINED_DESCRIPTION + "\n\n" +
            "<h4>Usage example:</h4>" +
            "<pre>" +
            "java -jar picard.jar DownsampleSam \\<br />" +
            "      I=input.bam \\<br />" +
            "      O=downsampled.bam" +
            "</pre>" +
            "<hr />";
    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to downsample.")
    public File INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output, downsampled, SAM or BAM file to write.")
    public File OUTPUT;

    @Option(shortName="S", doc="The downsampling strategy to use. See usage for discussion.")
    public Strategy STRATEGY = Strategy.ConstantMemory;

    @Option(shortName = "R", doc = "Random seed to use if deterministic behavior is desired.  " +
            "Setting to null will cause multiple invocations to produce different results.")
    public Integer RANDOM_SEED = 1;

    @Option(shortName = "P", doc = "The probability of keeping any individual read, between 0 and 1.")
    public double PROBABILITY = 1;

    @Option(shortName = "A", doc = "The accuracy that the downsampler should try to achieve if the selected strategy supports it. " +
            "Note that accuracy is never guaranteed, but some strategies will attempt to provide accuracy within the requested bounds." +
            "Higher accuracy will generally require more memory.")
    public double ACCURACY = 0.0001;

    @Option(shortName = "N", doc = "Keep this many individual reads.", mutex = "PROBABILITY")
    public long NUM_READS = -1;

    private final Log log = Log.getInstance(DownsampleSam.class);

    public static void main(final String[] args) {
        new DownsampleSam().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        if (0 < NUM_READS) {
            // NB: may count secondary or supplementary records, so it may be an overestimate of the number of records to be kept later.
            long numRecords = 0;
            final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
            if (in.hasIndex()) {
                final AbstractBAMFileIndex index = (AbstractBAMFileIndex) in.indexing().getIndex();
                final int nRefs = index.getNumberOfReferences();
                for (int i = 0; i < nRefs; i++) {
                    final BAMIndexMetaData indexMetaData = index.getMetaData(i);
                    numRecords += indexMetaData.getAlignedRecordCount() + indexMetaData.getUnalignedRecordCount();
                }
                numRecords += index.getNoCoordinateCount();
            }
            else {
                final Iterable<SAMRecord> iterable = in::iterator;
                numRecords = StreamSupport.stream(iterable.spliterator(), false).filter(record -> !record.isSecondaryOrSupplementary()).count();
            }
            if (NUM_READS > numRecords) {
                log.warn("Running DownsampleSam with NUM_READS greater than the number of reads in the input file.  Defaulting to 1.0.");
                PROBABILITY = 1.0;
            }
            else {
                PROBABILITY = NUM_READS / (double) numRecords;
            }
        }

        // Warn the user if they are running with P=1; 0 <= P <= 1 is checked by the DownsamplingIteratorFactory
        if (PROBABILITY == 1) {
            log.warn("Running DownsampleSam with PROBABILITY=1! This will likely just recreate the input file.");
        }

        final Integer randomSeed = RANDOM_SEED == null ? 42: RANDOM_SEED;
        final SamReader in = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);
        final SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(in.getFileHeader(), true, OUTPUT);
        final ProgressLogger progress = new ProgressLogger(log, (int) 1e7, "Wrote");
        final DownsamplingIterator iterator = DownsamplingIteratorFactory.make(INPUT, STRATEGY, PROBABILITY, ACCURACY, randomSeed);

        while (iterator.hasNext()) {
            final SAMRecord rec = iterator.next();
            out.addAlignment(rec);
            progress.record(rec);
        }

        out.close();
        CloserUtil.close(in);
        final NumberFormat fmt = new DecimalFormat("0.00%");
        log.info("Finished downsampling.");
        log.info("Kept ", iterator.getAcceptedCount(), " out of ", iterator.getSeenCount(), " reads (", fmt.format(iterator.getAcceptedFraction()), ").");

        return 0;
    }
}
