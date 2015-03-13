package picard.analysis.oxidation;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.InsertSizeFilter;
import htsjdk.samtools.filter.NotPrimaryAlignmentFilter;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Metrics;
import picard.util.DbSnpBitSetUtil;
import picard.analysis.oxidation.OxidationMetrics.*;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Original (non-SinglePassSamProgram) version.
 */
@CommandLineProgramProperties(
        usage = CollectDnaOxidationMetricsOriginal.USAGE,
        usageShort = CollectDnaOxidationMetricsOriginal.USAGE,
        programGroup = Metrics.class
)
public class CollectDnaOxidationMetricsOriginal extends CommandLineProgram {
    static final String USAGE = "Collect metrics relating to various kinds of DNA oxidation damage.";

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "Input BAM file for analysis.")
    public File INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc = "Base path of output files to write.")
    public File OUTPUT;

    @Option(shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME,
            doc = "Reference sequence to which BAM is aligned.")
    public File REFERENCE_SEQUENCE;

    @Option(doc = "An optional list of intervals to restrict analysis to.",
            optional = true)
    public File INTERVALS;

    @Option(doc = "VCF format dbSNP file, used to exclude regions around known polymorphisms from analysis.",
            optional = true)
    public File DB_SNP;

    @Option(shortName = "Q",
            doc = "The minimum base quality score for a base to be included in analysis.")
    public int MINIMUM_QUALITY_SCORE = 20;

    @Option(shortName = "MQ",
            doc = "The minimum mapping quality score for a base to be included in analysis.")
    public int MINIMUM_MAPPING_QUALITY = 30;

    @Option(shortName = "MIN_INS",
            doc = "The minimum insert size for a read to be included in analysis. Set of 0 to allow unpaired reads.")
    public int MINIMUM_INSERT_SIZE = 60;

    @Option(shortName = "MAX_INS",
            doc = "The maximum insert size for a read to be included in analysis. Set of 0 to allow unpaired reads.")
    public int MAXIMUM_INSERT_SIZE = 600;

    @Option(doc = "When available, use original quality scores for filtering.")
    public boolean USE_OQ = true;

    @Option(doc = "The number of context bases to include on each side of the assayed base. Note that the size of the output " +
            "metrics will grow exponentially with this number, so take care when increasing it.")
    public int CONTEXT_SIZE = 1;

    @Option(doc = "If specified, only print results for these contexts in the detail metrics output. The summary metrics output " +
            "will still take all contexts into consideration.")
    public Set<String> CONTEXTS_TO_PRINT = new HashSet<String>();

    @Option(doc = "For debugging purposes: stop after visiting this many sites with at least 1X coverage.")
    public long STOP_AFTER = Long.MAX_VALUE;

    private final Log log = Log.getInstance(CollectDnaOxidationMetricsOriginal.class);
    private static final String UNKNOWN_LIBRARY = "UnknownLibrary";
    private static final String UNKNOWN_SAMPLE = "UnknownSample";

    // Stock main method
    public static void main(final String[] args) {
        new CollectDnaOxidationMetricsOriginal().instanceMainWithExit(args);
    }

    @Override
    protected String[] customCommandLineValidation() {
        final int size = 1 + 2 * CONTEXT_SIZE;
        final List<String> messages = new ArrayList<String>();

        for (final String cxt : CONTEXTS_TO_PRINT) {
            for (final char base : cxt.toCharArray()) {
                if (!SequenceUtil.isValidBase((byte) base)) {
                    messages.add("Context " + cxt + " contains invalid bases");
                }
            }
            if (cxt.length() != size) {
                messages.add("Context " + cxt + " is not " + size + " long as implied by CONTEXT_SIZE=" + CONTEXT_SIZE);
            }
        }

        if (CONTEXT_SIZE < 0) messages.add("CONTEXT_SIZE cannot be negative");

        return messages.isEmpty() ? null : messages.toArray(new String[messages.size()]);
    }

    @Override
    protected int doWork() {
        final File PRE_ADAPT_SUMMARY_OUT = new File(OUTPUT + ".pre_adapt_summary_metrics");
        final File PRE_ADAPT_DETAILS_OUT = new File(OUTPUT + ".pre_adapt_detail_metrics");
        final File BAIT_BIAS_SUMMARY_OUT = new File(OUTPUT + ".bait_bias_summary_metrics");
        final File BAIT_BIAS_DETAILS_OUT = new File(OUTPUT + ".bait_bias_detail_metrics");

        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(PRE_ADAPT_SUMMARY_OUT);
        IOUtil.assertFileIsWritable(PRE_ADAPT_DETAILS_OUT);
        IOUtil.assertFileIsWritable(BAIT_BIAS_SUMMARY_OUT);
        IOUtil.assertFileIsWritable(BAIT_BIAS_DETAILS_OUT);

        if (INTERVALS != null) IOUtil.assertFileIsReadable(INTERVALS);
        IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);

        final ReferenceSequenceFileWalker refWalker = new ReferenceSequenceFileWalker(REFERENCE_SEQUENCE);
        final SamReader in = SamReaderFactory.makeDefault().open(INPUT);

        final Set<String> samples = new HashSet<String>();
        final Set<String> libraries = new HashSet<String>();
        for (final SAMReadGroupRecord rec : in.getFileHeader().getReadGroups()) {
            samples.add(nvl(rec.getSample(), UNKNOWN_SAMPLE));
            libraries.add(nvl(rec.getLibrary(), UNKNOWN_LIBRARY));
        }
        final String sampleAlias = StringUtil.join(",", new ArrayList<String>(samples));
        final int twoSidedContextLength = 2 * CONTEXT_SIZE + 1;

        // Load up dbSNP if available
        log.info("Loading dbSNP File: " + DB_SNP);
        final DbSnpBitSetUtil dbSnp;
        if (DB_SNP != null) dbSnp = new DbSnpBitSetUtil(DB_SNP, in.getFileHeader().getSequenceDictionary());
        else dbSnp = null;

        // Make an iterator that will filter out funny looking things
        final SamLocusIterator iterator;
        if (INTERVALS != null) {
            final IntervalList intervals = IntervalList.fromFile(INTERVALS);
            iterator = new SamLocusIterator(in, intervals.uniqued(), false);
        } else {
            iterator = new SamLocusIterator(in);
        }
        iterator.setEmitUncoveredLoci(false);
        iterator.setMappingQualityScoreCutoff(MINIMUM_MAPPING_QUALITY);
        iterator.setSamFilters(Arrays.asList(
                new NotPrimaryAlignmentFilter(),
                new DuplicateReadFilter(),
                new InsertSizeFilter(MINIMUM_INSERT_SIZE, MAXIMUM_INSERT_SIZE)
        ));

        // Set up the ArtifactCounters - one per library
        final Map<String, ArtifactCounter> artifactCounters = new HashMap<String, ArtifactCounter>();
        for (final String library : libraries) {
            artifactCounters.put(library, new ArtifactCounter(sampleAlias, library, CONTEXT_SIZE));
        }

        // Iterate over reference loci
        log.info("Starting iteration.");
        long nextLogTime = 0;
        int sites = 0;

        for (final SamLocusIterator.LocusInfo info : iterator) {
            // Skip dbSNP sites
            final String chrom = info.getSequenceName();
            final int pos = info.getPosition();
            final int index = pos - 1;
            if (dbSnp != null && dbSnp.isDbSnpSite(chrom, pos)) continue;

            // Skip sites at the end of chromosomes
            final byte[] bases = refWalker.get(info.getSequenceIndex()).getBases();
            if (pos < twoSidedContextLength || pos > bases.length - twoSidedContextLength) continue;

            // Get the reference context
            final String context = StringUtil.bytesToString(bases, index - CONTEXT_SIZE, twoSidedContextLength).toUpperCase();
            if (context.contains("N")) continue;

            // Iterate over records at this locus
            for (final SamLocusIterator.RecordAndOffset rao : info.getRecordAndPositions()) {
                final SAMRecord rec = rao.getRecord();
                final String library = nvl(rec.getReadGroup().getLibrary(), UNKNOWN_LIBRARY);
                final char readBase = Character.toUpperCase((char) rao.getReadBase());

                if (!libraries.contains(library)) continue;
                if (!passesBqCutoff(rao)) continue;
                if (readBase == 'N') continue;

                // count the base!
                artifactCounters.get(library).countRecord(context, readBase, rec);
            }

            // See if we need to stop
            ++sites;
            if (sites % 100 == 0) {
                final long now = System.currentTimeMillis();
                if (now > nextLogTime) {
                    log.info("Visited " + sites + " sites of interest. Last site: " + chrom + ":" + pos);
                    nextLogTime = now + 60000;
                }
            }
            if (sites >= STOP_AFTER) break;
        }

        // Finish up and write metrics
        final MetricsFile<PreAdaptSummaryMetrics, Integer> preAdaptSummaryMetricsFile = getMetricsFile();
        final MetricsFile<PreAdaptDetailMetrics, Integer> preAdaptDetailMetricsFile = getMetricsFile();
        final MetricsFile<BaitBiasSummaryMetrics, Integer> baitBiasSummaryMetricsFile = getMetricsFile();
        final MetricsFile<BaitBiasDetailMetrics, Integer> baitBiasDetailMetricsFile = getMetricsFile();

        for (final ArtifactCounter counter : artifactCounters.values()) {
            counter.finish();
            // write summary metrics
            preAdaptSummaryMetricsFile.addAllMetrics(counter.getPreAdaptSummaryMetrics());
            baitBiasSummaryMetricsFile.addAllMetrics(counter.getBaitBiasSummaryMetrics());
            // write detail metrics, subsetting if desired
            for (final PreAdaptDetailMetrics padm : counter.getPreAdaptDetailMetrics()) {
                if (CONTEXTS_TO_PRINT.size() == 0 || CONTEXTS_TO_PRINT.contains(padm.CONTEXT)) {
                    preAdaptDetailMetricsFile.addMetric(padm);
                }
            }
            for (final BaitBiasDetailMetrics bbdm : counter.getBaitBiasDetailMetrics()) {
                if (CONTEXTS_TO_PRINT.size() == 0 || CONTEXTS_TO_PRINT.contains(bbdm.CONTEXT)) {
                    baitBiasDetailMetricsFile.addMetric(bbdm);
                }
            }
        }

        preAdaptDetailMetricsFile.write(PRE_ADAPT_DETAILS_OUT);
        preAdaptSummaryMetricsFile.write(PRE_ADAPT_SUMMARY_OUT);
        baitBiasDetailMetricsFile.write(BAIT_BIAS_DETAILS_OUT);
        baitBiasSummaryMetricsFile.write(BAIT_BIAS_SUMMARY_OUT);

        CloserUtil.close(in);
        return 0;
    }

    private boolean passesBqCutoff(final SamLocusIterator.RecordAndOffset rao) {
        final byte qual;
        if (USE_OQ) {
            final byte[] oqs = rao.getRecord().getOriginalBaseQualities();
            if (oqs != null) qual = oqs[rao.getOffset()];
            else qual = rao.getBaseQuality();
        } else {
            qual = rao.getBaseQuality();
        }

        return (qual >= MINIMUM_QUALITY_SCORE);
    }

    /** Mimic of Oracle's nvl() - returns the first value if not null, otherwise the second value. */
    private <T> T nvl(final T value1, final T value2) {
        if (value1 != null) return value1;
        else return value2;
    }
}
