package picard.analysis.oxidation;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.SequenceUtil;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Metrics;
import picard.analysis.CollectOxoGMetrics.*;
import picard.analysis.oxidation.OxidationMetrics.*;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

@CommandLineProgramProperties(
        usage = ConvertDnaOxidationToOxoG.USAGE,
        usageShort = ConvertDnaOxidationToOxoG.USAGE,
        programGroup = Metrics.class
)
public class ConvertDnaOxidationToOxoG extends CommandLineProgram {
    static final String USAGE = "Extract OxoG metrics format from generalized DNA oxidation metrics.";

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "Basename for input DNA oxidation metrics")
    public File INPUT_BASE;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc = "Basename for output OxoG metrics. Defaults to same basename as input metrics",
            optional = true)
    public File OUTPUT_BASE;

    // Stock main method
    public static void main(final String[] args) {
        new ConvertDnaOxidationToOxoG().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        if (OUTPUT_BASE == null) { OUTPUT_BASE = INPUT_BASE; }

        final File PRE_ADAPT_IN = new File(INPUT_BASE + ".pre_adapt_detail_metrics");
        final File BAIT_BIAS_IN = new File(INPUT_BASE + ".bait_bias_detail_metrics");
        final File OXOG_OUT = new File(OUTPUT_BASE + ".oxog_metrics");

        IOUtil.assertFileIsReadable(PRE_ADAPT_IN);
        IOUtil.assertFileIsReadable(BAIT_BIAS_IN);
        IOUtil.assertFileIsWritable(OXOG_OUT);

        final List<PreAdaptDetailMetrics> preAdaptDetailMetricsList = (List<PreAdaptDetailMetrics>) MetricsFile.readBeans(PRE_ADAPT_IN);
        final List<BaitBiasDetailMetrics> baitBiasDetailMetricsList = (List<BaitBiasDetailMetrics>) MetricsFile.readBeans(BAIT_BIAS_IN);

        // TODO should we validate that the two inputs match up as expected?

        /**
         * Determine output fields. Just copy these from the input for now.
         */
        final String oxogSampleAlias = preAdaptDetailMetricsList.get(0).SAMPLE_ALIAS;
        final Set<String> oxogLibraries = new HashSet<String>();
        final Set<String> oxogContexts = new HashSet<String>();
        for (final PreAdaptDetailMetrics padm : preAdaptDetailMetricsList) {
            oxogLibraries.add(padm.LIBRARY);
            // Remember that OxoG only reports on the 'C' contexts
            // Also, exclude the summary row for each base (which is represented as a context of length 1)
            if (padm.REF_BASE == 'C' && padm.CONTEXT.length() > 1) {
                oxogContexts.add(padm.CONTEXT);
            }
        }

        /**
         * Store the inputs in maps of {Library -> {Context, Metric}} for easy access.
         * Remember, we only care about two transitions - C>A and G>T! Thus, for each context we
         * will only store one metric.
         */
        final Map<String, Map<String, PreAdaptDetailMetrics>> preAdaptDetailMetricsMap = new HashMap<String, Map<String, PreAdaptDetailMetrics>>();
        final Map<String, Map<String, BaitBiasDetailMetrics>> baitBiasDetailMetricsMap = new HashMap<String, Map<String, BaitBiasDetailMetrics>>();
        for (final String library : oxogLibraries) {
            final Map<String, PreAdaptDetailMetrics> contextsToPadm = new HashMap<String, PreAdaptDetailMetrics>();
            final Map<String, BaitBiasDetailMetrics> contextsToBbdm = new HashMap<String, BaitBiasDetailMetrics>();
            preAdaptDetailMetricsMap.put(library, contextsToPadm);
            baitBiasDetailMetricsMap.put(library, contextsToBbdm);
        }
        for (final PreAdaptDetailMetrics padm : preAdaptDetailMetricsList) {
            if (isOxoG(padm)) {
                preAdaptDetailMetricsMap.get(padm.LIBRARY).put(padm.CONTEXT, padm);
            }
        }
        for (final BaitBiasDetailMetrics bbdm : baitBiasDetailMetricsList) {
            if (isOxoG(bbdm)) {
                baitBiasDetailMetricsMap.get(bbdm.LIBRARY).put(bbdm.CONTEXT, bbdm);
            }
        }

        /**
         * Create the OxoG metrics
         */
        final List<CpcgMetrics> oxogMetrics = new ArrayList<CpcgMetrics>();
        for (final String library : oxogLibraries) {
            for (final String context : oxogContexts) {
                final CpcgMetrics m = new CpcgMetrics();
                m.SAMPLE_ALIAS = oxogSampleAlias;
                m.LIBRARY = library;
                m.CONTEXT = context;
                m.TOTAL_SITES = 0; // not calculated in the input metrics

                /**
                 * Get the relevant input metrics. This is done in a somewhat confusing way:
                 *
                 * 1. For pre-adapter metrics: note that OxoG only reports 'C' contexts, even though the actual pre-adapter OxoG artifact
                 *    occurs when the reference-strand base is 'G'. This is because OxoG reverse-complements all the contexts for some reason.
                 *    Thus when we add an entry for 'ACA' in the output, we actually need to get that data from 'TGT' in the input.
                 *
                 * 2. For bait-bias metrics: for each context, we report two opposing error rates, C_REF and G_REF, because for this metric
                 *    the bias could really go in either direction (whereas for pre-adapter artifacts we only expect one direction: G>T, but
                 *    never C>A, on the original reference strand). C_REF corresponds to the actual context printed in that row, and G_REF
                 *    corresponds to its reverse complement. So for 'ACA' in the output, we need to take data from both 'ACA' and 'TGT' in the
                 *    input.
                 */

                final PreAdaptDetailMetrics padm = preAdaptDetailMetricsMap.get(library).get(SequenceUtil.reverseComplement(context));
                final BaitBiasDetailMetrics bbdmFwd = baitBiasDetailMetricsMap.get(library).get(context);
                final BaitBiasDetailMetrics bbdmRev = baitBiasDetailMetricsMap.get(library).get(SequenceUtil.reverseComplement(context));

                // extract fields from PreAdaptDetailMetrics
                m.TOTAL_BASES = padm.PRO_REF_BASES + padm.PRO_ALT_BASES + padm.CON_REF_BASES + padm.CON_ALT_BASES;
                m.REF_TOTAL_BASES = padm.PRO_REF_BASES + padm.CON_REF_BASES;
                m.REF_NONOXO_BASES = padm.CON_REF_BASES;
                m.REF_OXO_BASES = padm.PRO_REF_BASES;
                m.ALT_NONOXO_BASES = padm.CON_ALT_BASES;
                m.ALT_OXO_BASES = padm.PRO_ALT_BASES;
                m.OXIDATION_ERROR_RATE = padm.ERROR_RATE;
                m.OXIDATION_Q = padm.QSCORE;

                // extract fields from BaitBiasDetailMetrics
                m.C_REF_REF_BASES = bbdmFwd.FWD_CXT_REF_BASES;
                m.G_REF_REF_BASES = bbdmFwd.REV_CXT_REF_BASES;
                m.C_REF_ALT_BASES = bbdmFwd.FWD_CXT_ALT_BASES;
                m.G_REF_ALT_BASES = bbdmFwd.REV_CXT_ALT_BASES;
                m.C_REF_OXO_ERROR_RATE = bbdmFwd.ERROR_RATE;
                m.C_REF_OXO_Q = bbdmFwd.QSCORE;
                m.G_REF_OXO_ERROR_RATE = bbdmRev.ERROR_RATE;
                m.G_REF_OXO_Q = bbdmRev.QSCORE;

                // add it
                oxogMetrics.add(m);
            }
        }

        final MetricsFile<CpcgMetrics, Integer> outputFile = getMetricsFile();
        for (final CpcgMetrics m : oxogMetrics) {
            outputFile.addMetric(m);
        }

        outputFile.write(OXOG_OUT);
        return 0;
    }

    /**
     * Convenience method for determining whether a metric is related to the OxoG artifact.
     */
    private boolean isOxoG(final DnaOxidationMetrics m) {
        final boolean cToA = (m.REF_BASE == 'C') && (m.ALT_BASE == 'A');
        final boolean gToT = (m.REF_BASE == 'G') && (m.ALT_BASE == 'T');
        return cToA || gToT;
    }
}
