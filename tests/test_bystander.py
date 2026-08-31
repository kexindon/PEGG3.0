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


def test_frame_source_from_row():
    """Per-row frame resolution: the strand and CDS travel with the variant, so a
    single run() call can mix transcripts. A row with no usable annotation must
    resolve to nothing rather than to a guess."""
    cds = [[100, 199], [300, 399]]

    ok = {'transcript_strand': '-', 'start_end_cds': cds, 'cds_valid': True}
    fmap, bounds, strand = by.frame_source_from_row(ok)
    assert strand == '-'
    assert fmap == by.cds_frame_map(cds, '-')
    assert bounds == by.cds_boundaries(cds)

    #a cache must not leak one transcript's frame into another
    cache = {}
    by.frame_source_from_row(ok, cache=cache)
    other = {'transcript_strand': '+', 'start_end_cds': cds, 'cds_valid': True}
    f2, _, s2 = by.frame_source_from_row(other, cache=cache)
    assert s2 == '+'
    assert f2 == by.cds_frame_map(cds, '+')
    assert f2 != fmap, 'strand must change the frame map'
    assert len(cache) == 2

    #every flavour of missing annotation yields no frame, never a guess
    import math
    for bad in ({'transcript_strand': '-', 'start_end_cds': cds, 'cds_valid': False},
                {'transcript_strand': None, 'start_end_cds': cds},
                {'transcript_strand': '-', 'start_end_cds': None},
                {'transcript_strand': float('nan'), 'start_end_cds': cds},
                {'transcript_strand': '-', 'start_end_cds': float('nan')},
                {'transcript_strand': '-', 'start_end_cds': []},
                {}):
        assert by.frame_source_from_row(bad) == (None, None, None), bad
    print('  per-row frame source                OK')


def test_attach_cds():
    """attach_cds() puts the annotation on the table in the columns run() reads."""
    import pandas as pd
    df = pd.DataFrame({'Hugo_Symbol': ['TP53', 'TP53', 'NOPE']})
    gene_cds = {'TP53': {'strand': '-', 'cds': [[100, 199]], 'valid': True}}
    out = by.attach_cds(df, gene_cds)

    assert list(out['transcript_strand']) == ['-', '-', None]
    assert list(out['cds_valid']) == [True, True, False]
    assert out['start_end_cds'][0] == [[100, 199]]
    assert out['start_end_cds'][2] is None
    #the unannotated row must resolve to no frame, so it gets ordinary pegRNAs
    assert by.frame_source_from_row(out.iloc[2]) == (None, None, None)
    #and the annotated one must resolve
    assert by.frame_source_from_row(out.iloc[0])[2] == '-'
    #input is not mutated
    assert 'start_end_cds' not in df.columns
    print('  attach_cds                          OK')


def test_exon_confinement():
    """A bystander must fall in the same exon as the edit.

    An RTT is one contiguous stretch of genomic DNA, so it cannot reach across an
    intron. Splice-buffer distance does not catch this on its own: two positions
    can each sit well inside their own exon and still be separated by an intron.
    """
    exons = [[1000, 1029], [1530, 1559]]          # 30 nt each, 500 nt intron
    frame_map = by.cds_frame_map(exons, '+')
    bounds = by.cds_boundaries(exons)
    positions_all = by.cds_positions(exons, '+')

    #an RTT straddling the junction: 12 nt of exon 1, then 12 nt of exon 2
    rtt_positions = positions_all[18:30] + positions_all[30:42]
    RTT = 'ATGGCTAGCACCGGTATGCTAGCA'
    left = 10                                      # edit 2 nt from exon 1's end

    def run(exon_blocks):
        return by.silent_bystanders(
            RTT, left, 1, 1, '+', '+', 0,
            RTT_genomic_positions=rtt_positions,
            frame_map=frame_map, boundaries=bounds,
            window_nt=9, max_muts=1, max_candidates=None,
            splice_buffer=0,        # isolate confinement from the splice buffer
            exon_blocks=exon_blocks)

    edit_genomic = rtt_positions[left]
    lo, hi = next((s, e) for s, e in exons if s <= edit_genomic <= e)

    def cross_exon(options):
        return sum(1 for o in options for p in o['positions']
                   if not (lo <= rtt_positions[p] <= hi))

    without = run(None)
    with_ = run(exons)

    #the bug is real: without confinement some options land in the next exon
    assert cross_exon(without) > 0, 'test is not exercising the junction'
    assert cross_exon(with_) == 0, 'bystander placed outside the edit\'s exon'
    assert len(with_) > 0, 'confinement should not remove every option here'
    print('  exon confinement                    OK (%d/%d cross-exon removed)'
          % (cross_exon(without), len(without)))


def test_exon_confinement_at_edge():
    """An edit right at an exon edge yields no bystanders rather than reaching
    into the neighbouring exon."""
    exons = [[1000, 1005], [1530, 1559]]          # tiny first exon
    frame_map = by.cds_frame_map(exons, '+')
    positions_all = by.cds_positions(exons, '+')
    rtt_positions = positions_all[:6] + positions_all[6:18]
    RTT = 'ATGGCTAGCACCGGTATGCT'[:len(rtt_positions)]

    options = by.silent_bystanders(
        RTT, 5, 1, 1, '+', '+', 0,
        RTT_genomic_positions=rtt_positions,
        frame_map=frame_map, boundaries=by.cds_boundaries(exons),
        window_nt=9, max_muts=1, max_candidates=None, splice_buffer=0,
        exon_blocks=exons)

    lo, hi = 1000, 1005
    for o in options:
        for p in o['positions']:
            assert lo <= rtt_positions[p] <= hi, 'escaped the edit exon'
    print('  exon confinement at an edge         OK (%d options)' % len(options))


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
    test_frame_source_from_row()
    test_attach_cds()
    test_exon_confinement()
    test_exon_confinement_at_edge()

    print()
    print('all tests passed')
