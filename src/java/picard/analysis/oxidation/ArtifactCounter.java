package picard.analysis.oxidation;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.ListMap;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import picard.PicardException;
import picard.analysis.oxidation.OxidationMetrics.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Keeps track of artifact counts, and extracts metrics once accumulation is finished.
 */
public class ArtifactCounter {
    private static final char N = 'N';

    private final String sampleAlias;
    private final String library;

    private final Set<String> fullContexts;
    private final Map<String, String> leadingContextMap;
    private final Map<String, String> trailingContextMap;
    private final Map<String, String> zeroContextMap;

    private final ContextAccumulator fullContextAccumulator;
    private final ContextAccumulator halfContextAccumulator;
    private final ContextAccumulator zeroContextAccumulator;

    private final List<PreAdaptSummaryMetrics> pasmList;
    private final List<PreAdaptDetailMetrics> padmList;
    private final List<BaitBiasSummaryMetrics> bbsmList;
    private final List<BaitBiasDetailMetrics> bbdmList;

    public ArtifactCounter(final String sampleAlias, final String library, final int contextSize) {
        this.sampleAlias = sampleAlias;
        this.library = library;

        // define the contexts
        this.fullContexts = new HashSet<String>();
        for (final byte[] kmer : SequenceUtil.generateAllKmers(2 * contextSize + 1)) {
            this.fullContexts.add(StringUtil.bytesToString(kmer));
        }

        // the half contexts specify either leading or trailing bases. the zero context is just the center.
        // NB: we use N to represent a wildcard base, rather than an ambiguous base. It's assumed that all of the input
        // contexts are unambiguous, and that any actual N's in the data have been dealt with elsewhere.
        final String padding = StringUtil.repeatCharNTimes(N, contextSize);
        this.leadingContextMap = new HashMap<String, String>();
        this.trailingContextMap = new HashMap<String, String>();
        this.zeroContextMap = new HashMap<String, String>();
        for (final String cxt : this.fullContexts) {
            final String leading = cxt.substring(0, contextSize);
            final String trailing = cxt.substring(contextSize + 1, cxt.length());
            final char center = cxt.charAt(contextSize);
            this.leadingContextMap.put(cxt, leading + center + padding);
            this.trailingContextMap.put(cxt, padding + center + trailing);
            this.zeroContextMap.put(cxt, padding + center + padding);
        }

        // set up the accumulators
        final Set<String> halfContexts = new HashSet<String>();
        halfContexts.addAll(leadingContextMap.values());
        halfContexts.addAll(trailingContextMap.values());
        final Set<String> zeroContexts = new HashSet<String>();
        zeroContexts.addAll(zeroContextMap.values());

        this.fullContextAccumulator = new ContextAccumulator(fullContexts);
        this.halfContextAccumulator = new ContextAccumulator(halfContexts);
        this.zeroContextAccumulator = new ContextAccumulator(zeroContexts);

        // these will get populated in the final step
        pasmList = new ArrayList<PreAdaptSummaryMetrics>();
        padmList = new ArrayList<PreAdaptDetailMetrics>();
        bbsmList = new ArrayList<BaitBiasSummaryMetrics>();
        bbdmList = new ArrayList<BaitBiasDetailMetrics>();
    }

    /**
     * Add a record to all the accumulators.
     */
    public void countRecord(final String refContext, final char calledBase, final SAMRecord rec) {
        this.fullContextAccumulator.countRecord(refContext, calledBase, rec);
        this.halfContextAccumulator.countRecord(this.leadingContextMap.get(refContext), calledBase, rec);
        this.halfContextAccumulator.countRecord(this.trailingContextMap.get(refContext), calledBase, rec);
        this.zeroContextAccumulator.countRecord(this.zeroContextMap.get(refContext), calledBase, rec);
    }

    /**
     * Stop counting, tally things up, and extract metrics.
     */
    public void finish() {
        final ListMap<Transition, DetailPair> allDetailMetrics = getDetailMetrics();
        final Map<Transition, SummaryPair> allSummaryMetrics = getSummaryMetrics();

        for (final Transition txn : Transition.altValues()) {
            final SummaryPair sp = allSummaryMetrics.get(txn);
            final List<DetailPair> dps = allDetailMetrics.get(txn);
            pasmList.add(sp.pasm);
            bbsmList.add(sp.bbsm);
            for (final DetailPair dp : dps) {
                padmList.add(dp.padm);
                bbdmList.add(dp.bbdm);
            }
        }
    }

    public List<PreAdaptSummaryMetrics> getPreAdaptSummaryMetrics() { return pasmList; }
    public List<PreAdaptDetailMetrics> getPreAdaptDetailMetrics() { return padmList; }
    public List<BaitBiasSummaryMetrics> getBaitBiasSummaryMetrics() { return bbsmList; }
    public List<BaitBiasDetailMetrics> getBaitBiasDetailMetrics() { return bbdmList; }

    /**
     * Core method to compute summary metrics. For each transition, we report:
     * 1. the total Q-score across all contexts
     * 2. the worst full context and its Q-score
     * 3. the worst leading context and its Q-score
     * 4. the worst trailing context and its Q-score
     *
     */
    private Map<Transition, SummaryPair> getSummaryMetrics() {
        final Map<Transition, SummaryPair> summaryMetricsMap = new HashMap<Transition, SummaryPair>();

        // extract the detail metrics from each accumulator
        final ListMap<Transition, DetailPair> fullMetrics = this.fullContextAccumulator.calculateMetrics(sampleAlias, library);
        final ListMap<Transition, DetailPair> halfMetrics = this.halfContextAccumulator.calculateMetrics(sampleAlias, library);
        final ListMap<Transition, DetailPair> zeroMetrics = this.zeroContextAccumulator.calculateMetrics(sampleAlias, library);

        // compute the summary metrics - one row for each transition
        for (final Transition txn : Transition.altValues()) {
            final List<DetailPair> fullMetricsForTxn = fullMetrics.get(txn);
            final List<DetailPair> zeroMetricsForTxn = zeroMetrics.get(txn);
            if (zeroMetricsForTxn.size() != 1) {
                throw new PicardException("Should have exactly one context-free metric pair for transition: " + txn);
            }

            // we want to report on leading / trailing contexts separately
            final List<DetailPair> leadingMetricsForTxn = new ArrayList<DetailPair>();
            final List<DetailPair> trailingMetricsForTxn= new ArrayList<DetailPair>();
            for (final DetailPair metrics : halfMetrics.get(txn)) {
                // first make sure they're the same context
                if (!metrics.padm.CONTEXT.equals(metrics.bbdm.CONTEXT)) {
                    throw new PicardException("Input detail metrics are not matched up properly - contexts differ.");
                }
                final boolean isLeading = this.leadingContextMap.containsValue(metrics.padm.CONTEXT);
                final boolean isTrailing = this.trailingContextMap.containsValue(metrics.padm.CONTEXT);
                // if the original contextSize is 0, there's no difference between leading and trailing, so add it to both
                if (isLeading) leadingMetricsForTxn.add(metrics);
                if (isTrailing) trailingMetricsForTxn.add(metrics);
            }

            // get the worst cases
            final DetailPair totalMetric = zeroMetricsForTxn.get(0);
            final DetailPair worstFullMetric = getWorstMetrics(fullMetricsForTxn);
            final DetailPair worstLeadingMetric = getWorstMetrics(leadingMetricsForTxn);
            final DetailPair worstTrailingMetric = getWorstMetrics(trailingMetricsForTxn);

            // construct the actual summary metrics - a combination of all the data we've just extracted
            final PreAdaptSummaryMetrics pasm = new PreAdaptSummaryMetrics();
            final BaitBiasSummaryMetrics bbsm = new BaitBiasSummaryMetrics();

            pasm.SAMPLE_ALIAS = this.sampleAlias;
            pasm.LIBRARY = this.library;
            pasm.REF_BASE = txn.ref();
            pasm.ALT_BASE = txn.call();
            pasm.TOTAL_QSCORE = totalMetric.padm.QSCORE;
            pasm.WORST_CXT = worstFullMetric.padm.CONTEXT;
            pasm.WORST_CXT_QSCORE = worstFullMetric.padm.QSCORE;
            pasm.WORST_PRE_CXT = worstLeadingMetric.padm.CONTEXT;
            pasm.WORST_PRE_CXT_QSCORE = worstLeadingMetric.padm.QSCORE;
            pasm.WORST_POST_CXT = worstTrailingMetric.padm.CONTEXT;
            pasm.WORST_POST_CXT_QSCORE = worstTrailingMetric.padm.QSCORE;
            pasm.inferArtifactName();

            bbsm.SAMPLE_ALIAS = this.sampleAlias;
            bbsm.LIBRARY = this.library;
            bbsm.REF_BASE = txn.ref();
            bbsm.ALT_BASE = txn.call();
            bbsm.TOTAL_QSCORE = totalMetric.bbdm.QSCORE;
            bbsm.WORST_CXT = worstFullMetric.bbdm.CONTEXT;
            bbsm.WORST_CXT_QSCORE = worstFullMetric.bbdm.QSCORE;
            bbsm.WORST_PRE_CXT = worstLeadingMetric.bbdm.CONTEXT;
            bbsm.WORST_PRE_CXT_QSCORE = worstLeadingMetric.bbdm.QSCORE;
            bbsm.WORST_POST_CXT = worstTrailingMetric.bbdm.CONTEXT;
            bbsm.WORST_POST_CXT_QSCORE = worstTrailingMetric.bbdm.QSCORE;
            bbsm.inferArtifactName();

            // add the finalized metrics to the map
            summaryMetricsMap.put(txn, new SummaryPair(pasm, bbsm));
        }
        return summaryMetricsMap;
    }

    private ListMap<Transition, DetailPair> getDetailMetrics() {
        return this.fullContextAccumulator.calculateMetrics(this.sampleAlias, this.library);
    }

    /**
     * Given a list of detail metrics, get the worst pre-adapter metrics, and independently from that get the worst bait bias metrics
     * (in terms of Q-score).
     */
    private DetailPair getWorstMetrics(final List<DetailPair> metrics) {
        PreAdaptDetailMetrics worstPadm = null;
        BaitBiasDetailMetrics worstBbdm = null;
        for (final DetailPair m : metrics) {
            if (worstPadm == null || m.padm.QSCORE < worstPadm.QSCORE) worstPadm = m.padm;
            if (worstBbdm == null || m.bbdm.QSCORE < worstBbdm.QSCORE) worstBbdm = m.bbdm;
        }
        return new DetailPair(worstPadm, worstBbdm);
    }
}
