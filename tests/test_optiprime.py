"""
Tests for pegg.optiprime's input construction.

These do NOT need OptiPrime installed: they check the rows pegg hands to the
model, which is where the bug these guard against lived. Scoring itself is
exercised separately (and needs the OptiPrime environment).

The expected lengths here are computed by hand from the design geometry rather
than from the module's own helpers, for the same reason as in test_bystander.
"""

import sys
import os

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from pegg import optiprime as op


#--- helpers ---------

def make_row(rtt_len, ref_len, alt_len, proto30=None):
    """
    A single pegRNA row shaped like pegg's output, with enough context either
    side that _rows_from_pegdf can place the protospacer window.

    The protospacer window must appear exactly once in the context, so the
    flanks are deliberately non-repetitive.
    """
    #4 nt pad + 20 nt protospacer + 6 nt (PAM plus downstream) = 30
    proto = proto30 if proto30 is not None else ('ACAC'
                                                 + 'CGTACGTTCGATCGTAACGT'
                                                 + 'CGGTAC')
    assert proto30 is not None or len(proto) == 30, len(proto)

    #long enough that the 25+RTT_length truncation never runs out of context
    left = 'TTTTGGGGAAAACCCCTTTTGGGGAAAACCCC' * 3
    right = 'GGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAA' * 3

    ref = 'A' * ref_len
    alt = 'T' * alt_len

    wt = left + proto + ref + right
    alt_seq = left + proto + alt + right

    return {
        'Protospacer_30': proto,
        'RTT': 'G' * rtt_len,
        'RTT_length': rtt_len,
        'PBS': 'C' * 13,
        'REF': ref,
        'ALT': alt,
        'wt_w_context': wt,
        'alt_w_context': alt_seq,
        'PAM_strand': '+',
    }


#--- tests ---------

def test_row_tuple_shape():
    """Every emitted row carries the six fields score() puts in the DataFrame."""
    df = pd.DataFrame([make_row(20, 1, 1)])
    rows, errors = op._rows_from_pegdf(df)
    assert len(rows) == 1 and len(errors) == 1, 'one row in, one row out'
    r = rows[0]
    assert r is not None, 'substitution row should be usable: %s' % errors[0]
    assert len(r) == 6, 'row must have 6 fields, got %d' % len(r)
    spacer, rtt, pbs, unedited, edited, proto30 = r
    assert len(spacer) == 20, 'spacer must be 20 nt, got %d' % len(spacer)
    assert len(proto30) == 30, 'proto30 must be 30 nt, got %d' % len(proto30)
    print('  row tuple shape                     OK')


def test_proto30_is_passed_through():
    """
    proto30 must be pegg's own Protospacer_30, not something the model derives.

    optiprime's format_pe_df otherwise computes proto30 = full_unedited[:30] and
    asserts it is exactly 30 nt. full_unedited is truncated to 25+RTT_length and,
    on an insertion, is shorter than full_edited by the inserted length -- so a
    short RTT drops it under 30 and the assertion kills the worker process,
    taking every pegRNA in that shard with it.
    """
    proto = 'ACAC' + 'CGTACGTTCGATCGTAACGT' + 'CGGTAC'
    df = pd.DataFrame([make_row(20, 1, 1, proto30=proto)])
    rows, _ = op._rows_from_pegdf(df)
    assert rows[0][5] == proto, 'proto30 must be carried through verbatim'
    assert rows[0][0] == proto[op.PS20_OFFSET:op.PS20_OFFSET + 20], \
        'spacer must be the 20-mer inside proto30'
    print('  proto30 passed through verbatim     OK')


def test_insertion_short_rtt_survives():
    """
    The regression: an insertion with a short RTT.

    full_unedited = 25 + RTT_length - inserted_length. With RTT=5 and a 3 nt
    insertion that is 27 -- under the 30 nt the model asserts on. The row must
    still be emitted, and must carry a 30 nt proto30 so the assertion holds.
    """
    bad = 0
    for rtt_len, ins_len in [(5, 3), (5, 1), (10, 9), (6, 6)]:
        df = pd.DataFrame([make_row(rtt_len, 1, 1 + ins_len)])
        rows, errors = op._rows_from_pegdf(df)
        r = rows[0]
        if r is None:
            print('     RTT=%-3d ins=%-3d NOT EMITTED (%s)' % (rtt_len, ins_len, errors[0]))
            bad += 1
            continue
        unedited, proto30 = r[3], r[5]
        if len(proto30) != 30:
            print('     RTT=%-3d ins=%-3d proto30=%d nt' % (rtt_len, ins_len, len(proto30)))
            bad += 1
        #the condition that used to crash the worker
        if len(unedited) < 30:
            print('     RTT=%-3d ins=%-3d full_unedited=%d nt (would have crashed '
                  'the old derivation)' % (rtt_len, ins_len, len(unedited)))
    assert bad == 0, '%d insertion cases still broken' % bad
    print('  insertion + short RTT               OK')


def test_short_protospacer_is_rejected_not_fatal():
    """
    A Protospacer_30 that is not 30 nt (context ran out at a contig edge) must
    be reported as a per-row error, never passed on to break a whole shard.
    """
    df = pd.DataFrame([make_row(20, 1, 1, proto30='ACGT' * 6)])  # 24 nt
    rows, errors = op._rows_from_pegdf(df)
    assert rows[0] is None, 'a 24 nt protospacer must not be emitted'
    assert errors[0] is not None and '30' in errors[0], \
        'the error should say what was wrong, got %r' % errors[0]
    print('  short Protospacer_30 rejected       OK')


def test_valid_groups_are_checked():
    """
    An unrecognised group name is silently accepted by optiprime's rate loader
    (group_factor stays 0 and it only logs at INFO), producing a plausible but
    meaningless score. pegg has to reject it here instead.
    """
    assert op.DEFAULT_GROUP in op.VALID_GROUPS, 'the default must be valid'
    assert len(op.VALID_GROUPS) == 12, \
        'expected 12 groups, got %d' % len(op.VALID_GROUPS)
    df = pd.DataFrame([make_row(20, 1, 1)])
    try:
        op.score(df, group='BOGUS_CELL')
    except op.OptiPrimeError as e:
        assert 'BOGUS_CELL' in str(e), 'the error should name the bad group'
        print('  invalid group rejected              OK')
        return
    except Exception as e:
        #no OptiPrime installed: the group check must still come first
        raise AssertionError('group check did not run before %s: %s'
                             % (type(e).__name__, e))
    raise AssertionError('a bogus group name was accepted')


if __name__ == '__main__':
    print('input construction:')
    test_row_tuple_shape()
    test_proto30_is_passed_through()
    test_insertion_short_rtt_survives()
    test_short_protospacer_is_rejected_not_fatal()
    test_valid_groups_are_checked()

    print()
    print('all tests passed')
