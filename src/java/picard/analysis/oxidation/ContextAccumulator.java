package picard.analysis.oxidation;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.ListMap;
import htsjdk.samtools.util.SequenceUtil;
import picard.PicardException;
import picard.analysis.oxidation.OxidationMetrics.*;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * Keeps track of the AlignmentAccumulators for each artifact / context of interest.
 */
public class ContextAccumulator {
    private static final char[] BASES = {'A', 'C', 'G', 'T'};

    private final Map<Transition, Map<String, AlignmentAccumulator>> artifactMap;

    public ContextAccumulator(final Set<String> contexts) {
        this.artifactMap = new HashMap<Transition, Map<String, AlignmentAccumulator>>();
        for (final Transition txn : Transition.values()) {
            this.artifactMap.put(txn, new HashMap<String, AlignmentAccumulator>());
        }
        for (final String cxt : contexts) {
            final char refBase = getCentralBase(cxt);
            for (final char calledBase : BASES) {
                final Transition txn = Transition.transitionOf(refBase, calledBase);
                this.artifactMap.get(txn).put(cxt, new AlignmentAccumulator());
            }
        }
    }

    public void countRecord(final String refContext, final char calledBase, final SAMRecord rec) {
        final char refBase = getCentralBase(refContext);
        final Transition txn = Transition.transitionOf(refBase, calledBase);
        this.artifactMap.get(txn).get(refContext).countRecord(rec);
    }

    /**
     * Core method to compute detailed (i.e. context-by-context) metrics from this accumulator.
     */
    public ListMap<Transition, DetailPair> calculateMetrics(final String sampleAlias, final String library) {
        final ListMap<Transition, DetailPair> detailMetricsMap = new ListMap<Transition, DetailPair>();
        for (final Transition txnAlt : Transition.altValues()) {
            final Transition txnRef = txnAlt.matchingRef();
            for (final String cxt : this.artifactMap.get(txnAlt).keySet()) {
                // each combination of artifact + context represents a single metric row
                final PreAdaptDetailMetrics padm = new PreAdaptDetailMetrics();
                final BaitBiasDetailMetrics bbdm = new BaitBiasDetailMetrics();

                // populate basic fields
                padm.SAMPLE_ALIAS = sampleAlias;
                padm.LIBRARY = library;
                padm.CONTEXT = cxt;
                padm.REF_BASE = txnAlt.ref();
                padm.ALT_BASE = txnAlt.call();

                bbdm.SAMPLE_ALIAS = sampleAlias;
                bbdm.LIBRARY = library;
                bbdm.CONTEXT = cxt;
                bbdm.REF_BASE = txnAlt.ref();
                bbdm.ALT_BASE = txnAlt.call();

                // retrieve all the necessary alignment counters.
                final AlignmentAccumulator fwdRefAlignments = this.artifactMap.get(txnRef).get(cxt);
                final AlignmentAccumulator fwdAltAlignments = this.artifactMap.get(txnAlt).get(cxt);
                final AlignmentAccumulator revRefAlignments = this.artifactMap.get(txnRef.complement()).get(SequenceUtil.reverseComplement(cxt));
                final AlignmentAccumulator revAltAlignments = this.artifactMap.get(txnAlt.complement()).get(SequenceUtil.reverseComplement(cxt));

                // do the actual tallying, as explained in the metric definitions
                padm.PRO_REF_BASES = fwdRefAlignments.R1_POS + fwdRefAlignments.R2_NEG + revRefAlignments.R1_NEG + revRefAlignments.R2_POS;
                padm.PRO_ALT_BASES = fwdAltAlignments.R1_POS + fwdAltAlignments.R2_NEG + revAltAlignments.R1_NEG + revAltAlignments.R2_POS;
                padm.CON_REF_BASES = fwdRefAlignments.R1_NEG + fwdRefAlignments.R2_POS + revRefAlignments.R1_POS + revRefAlignments.R2_NEG;
                padm.CON_ALT_BASES = fwdAltAlignments.R1_NEG + fwdAltAlignments.R2_POS + revAltAlignments.R1_POS + revAltAlignments.R2_NEG;

                bbdm.FWD_CXT_REF_BASES = fwdRefAlignments.R1_POS + fwdRefAlignments.R1_NEG + fwdRefAlignments.R2_POS + fwdRefAlignments.R2_NEG;
                bbdm.FWD_CXT_ALT_BASES = fwdAltAlignments.R1_POS + fwdAltAlignments.R1_NEG + fwdAltAlignments.R2_POS + fwdAltAlignments.R2_NEG;
                bbdm.REV_CXT_REF_BASES = revRefAlignments.R1_POS + revRefAlignments.R1_NEG + revRefAlignments.R2_POS + revRefAlignments.R2_NEG;
                bbdm.REV_CXT_ALT_BASES = revAltAlignments.R1_POS + revAltAlignments.R1_NEG + revAltAlignments.R2_POS + revAltAlignments.R2_NEG;

                // calculate error rates + Q-scores
                padm.calculateDerivedStatistics();
                bbdm.calculateDerivedStatistics();

                // add the finalized metrics to the map
                detailMetricsMap.add(txnAlt, new DetailPair(padm, bbdm));
            }
        }
        return detailMetricsMap;
    }

    private char getCentralBase(final String context) {
        if (context.length() % 2 == 0) throw new PicardException("Contexts cannot have an even number of bases: " + context);
        else return context.charAt(context.length() / 2);
    }

    /**
     * Little class for breaking down alignments by read1/read2 and positive/negative strand.
     */
    private static class AlignmentAccumulator {
        private long R1_POS = 0;
        private long R1_NEG = 0;
        private long R2_POS = 0;
        private long R2_NEG = 0;

        private void countRecord(final SAMRecord rec) {
            final boolean isNegativeStrand = rec.getReadNegativeStrandFlag();
            final boolean isReadTwo = rec.getReadPairedFlag() && rec.getSecondOfPairFlag();
            if (isReadTwo) {
                if (isNegativeStrand) this.R2_NEG++;
                else this.R2_POS++;
            } else {
                if (isNegativeStrand) this.R1_NEG++;
                else this.R1_POS++;
            }
        }
    }
}
