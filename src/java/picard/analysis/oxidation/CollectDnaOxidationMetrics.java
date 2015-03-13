package picard.analysis.oxidation;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.AggregateFilter;
import htsjdk.samtools.filter.AlignedFilter;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.FailsVendorReadQualityFilter;
import htsjdk.samtools.filter.InsertSizeFilter;
import htsjdk.samtools.filter.MappingQualityFilter;
import htsjdk.samtools.filter.NotPrimaryAlignmentFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.IntervalListReferenceSequenceMask;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.StringUtil;
import picard.analysis.SinglePassSamProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
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
 * Quantify substitution errors caused by mismatched base pairings during various
 * stages of sample / library prep. These are most commonly caused by oxidation
 * (e.g. the 8-oxo-G error mode), but could have other causes as well.
 *
 * We measure two distinct error types - artifacts that are introduced before
 * the addition of the read1/read2 adapters ("pre adapt") and those that are
 * introduced after target selection ("bait bias"). For each of these, we provide
 * summary metrics as well as detail metrics broken down by reference context
 * (the ref bases surrounding the substitution event).
 *
 * For a deeper explanation, see Costello et al. 2013:
 * http://www.ncbi.nlm.nih.gov/pubmed/23303777
 *
 * @author mattsooknah
 *
 */
@CommandLineProgramProperties(
        usage = CollectDnaOxidationMetrics.USAGE,
        usageShort = CollectDnaOxidationMetrics.USAGE,
        programGroup = Metrics.class
)
public class CollectDnaOxidationMetrics extends SinglePassSamProgram {
    static final String USAGE = "Collect metrics relating to various kinds of DNA oxidation damage.";

    @Option(doc = "An optional list of intervals to restrict analysis to.", optional = true)
    public File INTERVALS;

    @Option(doc = "VCF format dbSNP file, used to exclude regions around known polymorphisms from analysis.", optional = true)
    public File DB_SNP;

    @Option(shortName = "Q", doc = "The minimum base quality score for a base to be included in analysis.")
    public int MINIMUM_QUALITY_SCORE = 20;

    @Option(shortName = "MQ", doc = "The minimum mapping quality score for a base to be included in analysis.")
    public int MINIMUM_MAPPING_QUALITY = 30;

    @Option(shortName = "MIN_INS", doc = "The minimum insert size for a read to be included in analysis. Set of 0 to allow unpaired reads.")
    public int MINIMUM_INSERT_SIZE = 60;

    @Option(shortName = "MAX_INS", doc = "The maximum insert size for a read to be included in analysis. Set of 0 to allow unpaired reads.")
    public int MAXIMUM_INSERT_SIZE = 600;

    @Option(doc = "When available, use original quality scores for filtering.")
    public boolean USE_OQ = true;

    @Option(doc = "The number of context bases to include on each side of the assayed base.")
    public int CONTEXT_SIZE = 1;

    @Option(doc = "If specified, only print results for these contexts in the detail metrics output. " +
                  "However, the summary metrics output will still take all contexts into consideration.")
    public Set<String> CONTEXTS_TO_PRINT = new HashSet<String>();

    private static final String UNKNOWN_LIBRARY = "UnknownLibrary";
    private static final String UNKNOWN_SAMPLE = "UnknownSample";

    private final Log log = Log.getInstance(CollectDnaOxidationMetrics.class);

    private File preAdaptSummaryOut;
    private File preAdaptDetailsOut;
    private File baitBiasSummaryOut;
    private File baitBiasDetailsOut;

    private IntervalListReferenceSequenceMask intervalMask;
    private DbSnpBitSetUtil dbSnpMask;
    private SamRecordFilter recordFilter;

    private Set<String> samples = new HashSet<String>();
    private Set<String> libraries = new HashSet<String>();

    final Map<String, ArtifactCounter> artifactCounters = new HashMap<String, ArtifactCounter>();

    public static void main(final String[] args) {
        new CollectDnaOxidationMetrics().instanceMainWithExit(args);
    }

    @Override
    protected String[] customCommandLineValidation() {
        final List<String> messages = new ArrayList<String>();
        if (CONTEXT_SIZE < 0) messages.add("CONTEXT_SIZE cannot be negative");
        return messages.isEmpty() ? null : messages.toArray(new String[messages.size()]);
    }

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        preAdaptSummaryOut = new File(OUTPUT + ".pre_adapt_summary_metrics");
        preAdaptDetailsOut = new File(OUTPUT + ".pre_adapt_detail_metrics");
        baitBiasSummaryOut = new File(OUTPUT + ".bait_bias_summary_metrics");
        baitBiasDetailsOut = new File(OUTPUT + ".bait_bias_detail_metrics");

        IOUtil.assertFileIsWritable(preAdaptSummaryOut);
        IOUtil.assertFileIsWritable(preAdaptDetailsOut);
        IOUtil.assertFileIsWritable(baitBiasSummaryOut);
        IOUtil.assertFileIsWritable(baitBiasDetailsOut);

        for (final SAMReadGroupRecord rec : header.getReadGroups()) {
            samples.add(nvl(rec.getSample(), UNKNOWN_SAMPLE));
            libraries.add(nvl(rec.getLibrary(), UNKNOWN_LIBRARY));
        }

        if (INTERVALS != null) {
            IOUtil.assertFileIsReadable(INTERVALS);
            intervalMask = new IntervalListReferenceSequenceMask(IntervalList.fromFile(INTERVALS).uniqued());
        }

        if (DB_SNP != null) {
            IOUtil.assertFileIsReadable(DB_SNP);
            dbSnpMask = new DbSnpBitSetUtil(DB_SNP, header.getSequenceDictionary());
        }

        // set record-level filters
        recordFilter = new AggregateFilter(Arrays.asList(
                new FailsVendorReadQualityFilter(),
                new NotPrimaryAlignmentFilter(),
                new DuplicateReadFilter(),
                new AlignedFilter(true), // discard unmapped reads
                new MappingQualityFilter(MINIMUM_MAPPING_QUALITY),
                new InsertSizeFilter(MINIMUM_INSERT_SIZE, MAXIMUM_INSERT_SIZE)
        ));

        // set up the artifact counters
        final String sampleAlias = StringUtil.join(",", new ArrayList<String>(samples));
        for (final String library : libraries) {
            artifactCounters.put(library, new ArtifactCounter(sampleAlias, library, CONTEXT_SIZE));
        }
    }

    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence ref) {
        // see if the whole read should be skipped
        if (recordFilter.filterOut(rec)) return;

        final String library = nvl(rec.getReadGroup().getLibrary(), UNKNOWN_LIBRARY);
        if (!libraries.contains(library)) return;

        // iterate over aligned positions
        for (final AlignmentBlock block : rec.getAlignmentBlocks()) {
            for (int offset = 0; offset < block.getLength(); offset++) {
                // remember, these are 1-based!
                final int readPos = block.getReadStart() + offset;
                final int refPos = block.getReferenceStart() + offset;

                /**
                 * Skip regions outside of intervals.
                 *
                 * NB: IntervalListReferenceSequenceMask.get() has side-effects which assume
                 * that successive ReferenceSequence's passed to this method will be in-order
                 * (e.g. it will break if you call acceptRead() with chr1, then chr2, then chr1
                 * again). So this only works if the underlying iteration is coordinate-sorted.
                 */
                if (intervalMask != null && !intervalMask.get(ref.getContigIndex(), refPos)) continue;

                // skip dbSNP sites
                if (dbSnpMask != null && dbSnpMask.isDbSnpSite(ref.getName(), refPos)) continue;

                // skip contexts with N bases, and those that clip the end of the reference
                final String context = getContextOrNull(refPos, ref);
                if (context == null || context.contains("N")) continue;

                // skip low BQ sites
                if (failsBaseQualityCutoff(readPos, rec)) continue;

                // skip N bases in read
                final char readBase = Character.toUpperCase((char) rec.getReadBases()[readPos - 1]);
                if (readBase == 'N') continue;

                // count the base!
                artifactCounters.get(library).countRecord(context, readBase, rec);
            }
        }
    }

    @Override
    protected void finish() {
        final MetricsFile<PreAdaptSummaryMetrics, Integer> preAdaptSummaryMetricsFile = getMetricsFile();
        final MetricsFile<PreAdaptDetailMetrics, Integer> preAdaptDetailMetricsFile = getMetricsFile();
        final MetricsFile<BaitBiasSummaryMetrics, Integer> baitBiasSummaryMetricsFile = getMetricsFile();
        final MetricsFile<BaitBiasDetailMetrics, Integer> baitBiasDetailMetricsFile = getMetricsFile();

        for (final ArtifactCounter counter : artifactCounters.values()) {
            // build metrics
            counter.finish();

            // write metrics
            preAdaptSummaryMetricsFile.addAllMetrics(counter.getPreAdaptSummaryMetrics());
            baitBiasSummaryMetricsFile.addAllMetrics(counter.getBaitBiasSummaryMetrics());

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

        preAdaptDetailMetricsFile.write(preAdaptDetailsOut);
        preAdaptSummaryMetricsFile.write(preAdaptSummaryOut);
        baitBiasDetailMetricsFile.write(baitBiasDetailsOut);
        baitBiasSummaryMetricsFile.write(baitBiasSummaryOut);
    }

    @Override
    protected boolean usesNoRefReads() { return false; }

    /**
     * Get the sequence context string (upper-case) at the specified position on the given reference.
     * Return null if the context is clipped by the end of the sequence.
     */
    private String getContextOrNull(final int oneIndexedPos, final ReferenceSequence ref) {
        final int contextStartIndex = oneIndexedPos - CONTEXT_SIZE - 1;
        final int contextFullLength = 2 * CONTEXT_SIZE + 1;
        try {
            return StringUtil.bytesToString(ref.getBases(), contextStartIndex, contextFullLength).toUpperCase();
        } catch (final IndexOutOfBoundsException e) {
            // catching the exception is perhaps sloppier than checking the bounds manually, but is less susceptible to bugs...
            return null;
        }
    }

    /**
     * Check if this read base fails the base quality cutoff.
     */
    private boolean failsBaseQualityCutoff(final int oneIndexedPos, final SAMRecord rec) {
        final byte qual;
        if (USE_OQ && rec.getOriginalBaseQualities() != null) {
            qual = rec.getOriginalBaseQualities()[oneIndexedPos - 1];
        } else {
            qual = rec.getBaseQualities()[oneIndexedPos - 1];
        }
        return (qual < MINIMUM_QUALITY_SCORE);
    }

    /** Mimic of Oracle's nvl() - returns the first value if not null, otherwise the second value. */
    private <T> T nvl(final T value1, final T value2) {
        if (value1 != null) return value1;
        else return value2;
    }
}
