"""
Tests for pegg.bystander.

The reading frames used here are counted by hand and hard-coded, rather than
derived with the module's own helpers. Deriving them would test the code against
its own assumptions, which is exactly how the frame-anchor bug survived its first
round of "validation".
"""

import sys
import os

import Bio.Seq

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from pegg import bystander as by


#--- helpers ---------

def translate_span(seq, lo, hi):
    """Translates seq[lo:hi], which must be codon-aligned and a multiple of 3."""
    span = seq[lo:hi]
    assert len(span) % 3 == 0, 'span is not a whole number of codons'
    return str(Bio.Seq.Seq(span).translate())


def transcript_view(RTT_fwd, transcript_strand, PAM_strand):
    """The RTT as the ribosome would read it."""
    if transcript_strand == PAM_strand:
        return RTT_fwd
    return str(Bio.Seq.Seq(RTT_fwd).reverse_complement())


def assert_option_valid(opt, RTT_fwd, left_RTT_len, ref_len, alt_len,
                        transcript_strand, PAM_strand, phase_of_work_start,
                        window_nt, max_muts):
    """
    Every invariant a returned option has to satisfy, checked independently of
    the module's internals. phase_of_work_start is the hand-counted codon phase
    of the first base of the transcript-orientation view.
    """
    L = len(RTT_fwd)

    #--- length is preserved
    assert len(opt['RTT']) == L, 'RTT length changed'

    #--- no more changes than asked for
    assert opt['n_muts'] <= max_muts, 'too many mutations'

    #--- the reported positions are the positions that actually differ
    actual = [i for i in range(L) if RTT_fwd[i].upper() != opt['RTT'][i].upper()]
    assert actual == opt['positions'], \
        'reported positions %s != actual %s' % (opt['positions'], actual)
    assert opt['n_muts'] == len(actual), 'n_muts disagrees with positions'

    #--- the intended edit is untouched
    if alt_len == 0:
        #a pure deletion has no alternate bases; the seam either side is protected
        forbidden = {left_RTT_len - 1, left_RTT_len}
    else:
        forbidden = set(range(left_RTT_len, left_RTT_len + alt_len))
    assert not (set(opt['positions']) & forbidden), \
        'option altered the intended edit at %s' % sorted(set(opt['positions']) & forbidden)

    #--- every change is synonymous, read in the transcript's own frame
    old_work = transcript_view(RTT_fwd, transcript_strand, PAM_strand)
    new_work = transcript_view(opt['RTT'], transcript_strand, PAM_strand)

    lo = (-phase_of_work_start) % 3
    hi = lo + 3 * ((L - lo) // 3)
    assert translate_span(old_work, lo, hi) == translate_span(new_work, lo, hi), \
        'option is NOT synonymous'


#--- fixture ---------

#A 23 nt RTT. Read on the + strand starting at offset 0 the codons are
#    ATG CCC GGG AAA TTT CCC GGT  ->  M P G K F P G
#which gives several degenerate codons (Pro, Gly, Lys, Phe) to work with.
RTT = 'ATGCCCGGGAAATTTCCCGGTAC'

#When the transcript runs on the opposite strand the ribosome reads
#    reverse_complement(RTT) = GTACCGGGAAATTTCCCGGGCAT
#Counted by hand: the RTT occupies p0..p22 on the PAM strand, and the transcript
#reads p22 -> p0. Taking frame_of_RTT_start = 1 (i.e. phase(p0) = 1), the phase
#decreases by one per base as we walk from p0 towards p22, so
#    phase(p22) = (1 - 22) mod 3 = 0
#and the reverse complement is codon-aligned from its very first base:
#    GTA CCG GGA AAT TTC CCG GGC | AT  ->  V P G N F P G
REV_FRAME_OF_RTT_START = 1
REV_PHASE_OF_WORK_START = 0

#On the same strand the anchor passes through untouched.
FWD_FRAME_OF_RTT_START = 0
FWD_PHASE_OF_WORK_START = 0

EDIT_OFFSET = 9      #the edit sits at RTT[9:], inside the AAA (Lys) codon
WINDOW = 5
MAX_MUTS = 2


#--- the six cases ---------

CASES = [
    #(label, transcript_strand, PAM_strand, ref_len, alt_len)
    ('fwd  substitution',          '+', '+', 1, 1),
    ('fwd  in-frame indel (+3)',   '+', '+', 0, 3),
    ('fwd  frameshift indel (+1)', '+', '+', 0, 1),
    ('rev  substitution',          '-', '+', 1, 1),
    ('rev  in-frame indel (-3)',   '-', '+', 3, 0),
    ('rev  frameshift indel (+2)', '-', '+', 0, 2),
]


def run_case(label, ts, ps, ref_len, alt_len, verbose=True):
    """Runs one case and asserts every invariant on every option returned."""
    if ts == ps:
        frame_arg = FWD_FRAME_OF_RTT_START
        work_phase = FWD_PHASE_OF_WORK_START
    else:
        frame_arg = REV_FRAME_OF_RTT_START
        work_phase = REV_PHASE_OF_WORK_START

    opts = by.silent_bystanders(
        RTT, EDIT_OFFSET, ref_len, alt_len, ts, ps, frame_arg,
        window_nt=WINDOW, max_muts=MAX_MUTS, max_candidates=None)

    for opt in opts:
        assert_option_valid(opt, RTT, EDIT_OFFSET, ref_len, alt_len,
                            ts, ps, work_phase, WINDOW, MAX_MUTS)

        #every change must sit strictly inside the window around the edit
        if alt_len == 0:
            lo, hi = EDIT_OFFSET - 1 - WINDOW, EDIT_OFFSET + 1 + WINDOW
        else:
            lo, hi = EDIT_OFFSET - WINDOW, EDIT_OFFSET + alt_len + WINDOW
        for p in opt['positions']:
            assert lo <= p < hi, \
                '%s: change at %d outside window [%d,%d)' % (label, p, lo, hi)

    #a frameshifting indel must not place anything downstream of the edit,
    #because the reading frame there no longer exists
    if (alt_len - ref_len) % 3 != 0:
        for opt in opts:
            for p in opt['positions']:
                if ts == ps:
                    assert p < EDIT_OFFSET, \
                        '%s: change at %d is downstream of a frameshift' % (label, p)
                else:
                    assert p >= EDIT_OFFSET, \
                        '%s: change at %d is downstream of a frameshift' % (label, p)

    if verbose:
        print('  %-28s options=%3d  OK' % (label, len(opts)))

    return opts


#--- regression tests for bugs that have actually occurred ---------

def test_reverse_frame_anchor():
    """
    The anchor must be translated when moving to the reverse complement.
    Using the untranslated anchor happens to be correct only when
    (L - 1) % 3 == 0, which is how the bug survived its first review.
    """
    assert by.reverse_frame_anchor(1, 23) == 0
    assert by.reverse_frame_anchor(0, 23) == 2
    assert by.reverse_frame_anchor(2, 23) == 1

    #the coincidence: unchanged exactly when (L-1) % 3 == 0
    for L in range(4, 40):
        same = [by.reverse_frame_anchor(f, L) == f for f in (0, 1, 2)]
        assert all(same) == ((L - 1) % 3 == 0), 'L=%d' % L
    print('  reverse_frame_anchor                OK')


def test_deletion_seam_protected():
    """
    A pure deletion has alt_len == 0, so the naive protected range is empty and
    the bases either side of the deletion seam were free to be changed.
    """
    opts = by.silent_bystanders(RTT, EDIT_OFFSET, 3, 0, '+', '+',
                                FWD_FRAME_OF_RTT_START,
                                window_nt=WINDOW, max_muts=MAX_MUTS,
                                max_candidates=None)
    for opt in opts:
        assert EDIT_OFFSET - 1 not in opt['positions']
        assert EDIT_OFFSET not in opt['positions']
    print('  deletion seam protected             OK (%d options)' % len(opts))


def test_wrong_anchor_is_not_caught_internally():
    """
    Documents a limitation worth knowing about: the internal translation check
    cannot catch a wrong frame anchor.

    The check re-translates in the frame the function was given. A wrong anchor
    is self-consistent -- the codons enumerated in the wrong frame really are
    synonymous in that wrong frame -- so the output passes the internal check
    while being non-synonymous in the true frame.

    The consequence is that reverse_frame_anchor() has no runtime safety net and
    is guarded only by tests. This test pins the failure mode so that nobody
    mistakes the internal check for protection it does not provide.
    """
    original = by.reverse_frame_anchor
    by.reverse_frame_anchor = lambda f, L: f      # the original bug
    try:
        opts = by.silent_bystanders(RTT, EDIT_OFFSET, 1, 1, '-', '+',
                                    REV_FRAME_OF_RTT_START,
                                    window_nt=WINDOW, max_muts=MAX_MUTS,
                                    max_candidates=None)
        non_synonymous = 0
        for opt in opts:
            try:
                assert_option_valid(opt, RTT, EDIT_OFFSET, 1, 1, '-', '+',
                                    REV_PHASE_OF_WORK_START, WINDOW, MAX_MUTS)
            except AssertionError:
                non_synonymous += 1
    finally:
        by.reverse_frame_anchor = original

    assert non_synonymous > 0, \
        'expected the broken anchor to produce non-synonymous output'
    print('  wrong anchor NOT caught internally  OK (%d/%d bad, tests are the guard)'
          % (non_synonymous, len(opts)))


def test_no_duplicate_rtts():
    for label, ts, ps, ref_len, alt_len in CASES:
        opts = run_case(label, ts, ps, ref_len, alt_len, verbose=False)
        seqs = [o['RTT'] for o in opts]
        assert len(seqs) == len(set(seqs)), '%s: duplicate RTTs' % label
    print('  no duplicate RTTs                   OK')


def test_max_candidates_caps():
    opts = by.silent_bystanders(RTT, EDIT_OFFSET, 1, 1, '+', '+',
                                FWD_FRAME_OF_RTT_START,
                                window_nt=WINDOW, max_muts=MAX_MUTS,
                                max_candidates=3)
    assert len(opts) <= 3
    print('  max_candidates caps                 OK')


def test_splice_and_cds_respected():
    """Positions outside the CDS or near an exon boundary are never altered."""
    cds = [[1000, 1008], [2000, 2020]]
    fm = by.cds_frame_map(cds, '+')
    bd = by.cds_boundaries(cds)

    #RTT laid across the junction: 9 exonic bases, then the next exon
    gpos = list(range(1000, 1009)) + list(range(2000, 2014))
    assert len(gpos) == len(RTT)

    opts = by.silent_bystanders(RTT, EDIT_OFFSET, 1, 1, '+', '+',
                                FWD_FRAME_OF_RTT_START,
                                RTT_genomic_positions=gpos, frame_map=fm,
                                boundaries=bd, window_nt=WINDOW,
                                max_muts=MAX_MUTS, max_candidates=None,
                                splice_buffer=3)
    for opt in opts:
        for p in opt['positions']:
            assert by.position_is_safe(gpos[p], fm, bd, 3), \
                'altered unsafe genomic position %d' % gpos[p]
    print('  CDS / splice buffer respected       OK (%d options)' % len(opts))


def test_codon_tables():
    import Bio.Seq as _S
    assert len(by.CODON_TO_AA) == 64
    for codon, aa in by.CODON_TO_AA.items():
        assert str(_S.Seq(codon).translate()) == aa
    assert by.synonymous_codons('ATG') == []
    assert by.synonymous_codons('TGG') == []
    assert by.translate_codon('NNN') is None
    print('  codon tables                        OK')


if __name__ == '__main__':
    print('six input classes:')
    for label, ts, ps, ref_len, alt_len in CASES:
        run_case(label, ts, ps, ref_len, alt_len)

    print()
    print('regression tests:')
    test_reverse_frame_anchor()
    test_deletion_seam_protected()
    test_wrong_anchor_is_not_caught_internally()
    test_no_duplicate_rtts()
    test_max_candidates_caps()
    test_splice_and_cds_respected()
    test_codon_tables()

    print()
    print('all tests passed')
