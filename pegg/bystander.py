"""
Silent bystander mutations for pegRNA design.

Adds synonymous ("silent") mutations into the RTT of a pegRNA, alongside the
intended edit. Silent changes near the edit help evade mismatch repair (MMR) and,
when they happen to fall in the PAM or seed, also reduce re-nicking of the edited
allele -- both of which can substantially increase prime editing efficiency.

This is conceptually similar to the silent-bystander addon of PRIDICT2.0
(https://github.com/uzh-dqbm-cmi/PRIDICT2): synonymous codon substitutions are
enumerated in a window around the edit. PRIDICT2.0 then ranks them with a deep
model; here they are handed back to PEGG's own pipeline so that they can be
generated at library scale.

Everything in this module is opt-in. Ordinary pegRNA design is unaffected by it
and needs no reading frame information; the requirements below apply only when
silent bystanders are requested, via prime.run(..., silent_bystander=True).

Reading frame
--------------
Every option returned is verified synonymous against an explicit reading frame,
which must be supplied by the caller. Two things make this less trivial than it
looks, and both are handled here:

(1) The RTT is stored in the PAM-strand orientation. When the transcript runs on
    the opposite strand, codons must be read on the reverse complement, and the
    frame anchor has to be translated along with the sequence -- the phase of
    work_seq[0] is not the phase of RTT_fwd[0]. See reverse_frame_anchor().

(2) Downstream of an insertion or deletion the reading frame is shifted, so
    "synonymous" is only well defined there if the net length change is a
    multiple of three. Codons that fall downstream of a frameshifting indel are
    excluded rather than guessed at.

Input formats
--------------
PEGG accepts three input formats. They do not all support silent bystanders
equally, because they do not all carry genomic coordinates:

    cBioPortal   Coordinates are available, so the reading frame is looked up
                 per position from a CDS annotation (start_end_cds +
                 transcript_strand). CDS membership and splice-site distance are
                 both enforced. Full support.

    WT_ALT       No coordinates. The reading frame must be declared with
    PrimeDesign  ORF_start (0, 1 or 2), and THE INPUT SEQUENCE MUST BE IN FRAME
                 -- entirely coding sequence, with the frame beginning at offset
                 ORF_start. CDS membership and splice sites cannot be checked in
                 this mode, so a sequence containing intronic or untranslated
                 bases would produce bystanders that are silent only in
                 appearance. Use cBioPortal input for anything spanning a splice
                 junction.

Again, none of this constrains ordinary pegRNA design: a sequence that is not in
frame, or not coding at all, is designed against exactly as before -- it simply
cannot carry silent bystanders. resolve_frame_source() enforces these
requirements rather than guessing, and only when the feature is switched on.
"""

from itertools import product

import Bio.Seq
from Bio.Data import CodonTable


#--- codon tables ---------

#Standard genetic code, taken from Biopython rather than hardcoded so that the
#table can never drift out of sync with the translations used elsewhere.
STANDARD_TABLE = CodonTable.unambiguous_dna_by_id[1]

#Codon -> amino acid. Biopython keeps stop codons out of forward_table, so they
#are folded back in here as '*' to match the convention used by Seq.translate().
CODON_TO_AA = dict(STANDARD_TABLE.forward_table)
for _stop in STANDARD_TABLE.stop_codons:
    CODON_TO_AA[_stop] = '*'

#Amino acid -> synonymous codons; inverted from CODON_TO_AA so the two agree.
CODON_MAP = {}
for _codon, _aa in CODON_TO_AA.items():
    CODON_MAP.setdefault(_aa, []).append(_codon)


def translate_codon(codon):
    """
    Translates a single codon to its amino acid. Returns None if the codon
    contains any non-ATCG character (e.g. an N in the reference genome).

    Parameters
    -----------
    codon
        *type = str*

        Three nucleotides.
    """
    return CODON_TO_AA.get(codon.upper())


def synonymous_codons(codon):
    """
    Returns the list of codons synonymous with the input codon, excluding the
    input codon itself. Returns an empty list for Met/Trp (no synonyms) and for
    codons containing non-ATCG characters.

    Parameters
    -----------
    codon
        *type = str*

        Three nucleotides.
    """
    aa = translate_codon(codon)
    if aa is None:
        return []
    return [i for i in CODON_MAP[aa] if i != codon.upper()]


#--- reading frame ---------

def cds_positions(start_end_cds, strand):
    """
    Flattens a transcript's CDS blocks into a list of genomic positions in
    transcription order (i.e. reversed for a - strand transcript).

    This is the single walk that the rest of the frame handling is built on:
    cds_frame_map() indexes it by position, and library.neutral_substitutions()
    chunks it into codons.

    Parameters
    -----------
    start_end_cds
        *type = list*

        A 2-d list containing the start/end locations of each region of the
        coding sequence (CDS) for the gene's selected transcript, ordered by
        position in the + strand orientation, 1-based and inclusive.
        e.g. [[7572930, 7573008], ...]

    strand
        *type = str*

        Strand that the transcript is on. Options are '+' or '-'.
    """
    assert strand in ['+', '-'], "strand must be '+' or '-'"

    positions = []
    for i in start_end_cds:
        for k in range(i[0], i[1] + 1):
            positions.append(k)

    #transcription order: - strand transcripts read from high to low coordinate
    if strand == '-':
        positions = positions[::-1]

    return positions


def cds_codons(start_end_cds, strand):
    """
    Groups a transcript's CDS positions into codons. Returns a list of
    [pos1, pos2, pos3] triples in transcription order. A trailing partial codon
    (i.e. a CDS whose length is not a multiple of three) is dropped.

    Parameters
    -----------
    start_end_cds
        *type = list*

        CDS blocks, 1-based inclusive, ordered in the + strand orientation.

    strand
        *type = str*

        Strand that the transcript is on. Options are '+' or '-'.
    """
    positions = cds_positions(start_end_cds, strand)

    codons = []
    for i in range(0, len(positions) - 2, 3):
        codons.append([positions[i], positions[i + 1], positions[i + 2]])

    return codons


def cds_frame_map(start_end_cds, strand):
    """
    Builds a map from genomic position to codon phase for a transcript's coding
    sequence. Returns a dict of {genomic_position: (codon_index, phase)} where
    codon_index counts codons from the start codon (0-based) and phase is the
    position within that codon (0, 1 or 2).

    Because the map is built by walking the CDS blocks in transcription order,
    codons split across an exon junction are handled correctly -- something a
    naive (position % 3) calculation cannot do.

    Parameters
    -----------
    start_end_cds
        *type = list*

        A 2-d list containing the start/end locations of each region of the
        coding sequence (CDS) for the gene's selected transcript, ordered by
        position in the + strand orientation. Same format as used by
        library.neutral_substitutions(). e.g. [[7572930, 7573008], ...]

    strand
        *type = str*

        Strand that the transcript is on. Options are '+' or '-'.
    """
    positions = cds_positions(start_end_cds, strand)

    frame_map = {}
    for idx, pos in enumerate(positions):
        frame_map[pos] = (idx // 3, idx % 3)

    return frame_map


def frame_at(genomic_position, frame_map):
    """
    Looks up the codon phase of a genomic position. Returns (codon_index, phase),
    or None if the position falls outside the coding sequence (e.g. in a UTR,
    intron, or a different gene).

    Parameters
    -----------
    genomic_position
        *type = int*

        Genomic coordinate to look up.

    frame_map
        *type = dict*

        Frame map generated by cds_frame_map().
    """
    return frame_map.get(genomic_position)



#--- splice site / CDS safety ---------

def cds_boundaries(start_end_cds):
    """
    Returns the set of genomic positions that sit at the edge of a CDS block
    (i.e. immediately adjacent to an exon/intron junction). Used to keep silent
    bystanders away from splice donor/acceptor sites.

    Parameters
    -----------
    start_end_cds
        *type = list*

        CDS blocks, 1-based inclusive, ordered in the + strand orientation.
        See cds_frame_map().
    """
    edges = set()
    for s, e in start_end_cds:
        edges.add(s)
        edges.add(e)
    return edges


def position_is_safe(genomic_position, frame_map, boundaries, splice_buffer=3):
    """
    True if a genomic position may be altered by a silent bystander mutation.
    A position is safe only if it lies within the CDS and is at least
    splice_buffer nt away from every exon boundary.

    Altering a base close to an exon boundary can destroy a splice donor or
    acceptor even when the change is synonymous at the protein level, which would
    cause exon skipping -- an effect that is invisible in the sequence but severe
    in the protein.

    Parameters
    -----------
    genomic_position
        *type = int*

        Genomic coordinate to test.

    frame_map
        *type = dict*

        Frame map generated by cds_frame_map().

    boundaries
        *type = set*

        Exon boundary positions generated by cds_boundaries().

    splice_buffer
        *type = int*

        Minimum distance to keep from an exon boundary. Default = 3.
    """
    if genomic_position not in frame_map:
        return False

    for b in boundaries:
        if abs(genomic_position - b) < splice_buffer:
            return False

    return True


#--- silent bystander search ---------

def reverse_frame_anchor(frame_of_RTT_start, RTT_length):
    """
    Translates a reading-frame anchor from the PAM strand to the reverse
    complement. Returns the codon phase of the first base of
    reverse_complement(RTT_fwd).

    Let the RTT occupy positions p_0 .. p_(L-1) read 5'->3' on the PAM strand,
    with phase(p_0) = frame_of_RTT_start. A transcript on the opposite strand is
    read in the other direction, p_(L-1) -> p_0, with the phase incrementing by
    one per base along that direction. The first base of the reverse complement
    is p_(L-1), so

        phase(p_(L-1)) = (frame_of_RTT_start - (L - 1)) % 3

    Note this equals frame_of_RTT_start only when (L - 1) % 3 == 0; using the
    untranslated anchor silently produces non-synonymous "silent" mutations for
    the other two thirds of RTT lengths.

    Parameters
    -----------
    frame_of_RTT_start
        *type = int*

        Codon phase (0, 1 or 2) of RTT_fwd[0], read in the transcript's
        orientation.

    RTT_length
        *type = int*

        Length of the RTT in nt.
    """
    return (frame_of_RTT_start - (RTT_length - 1)) % 3


def silent_bystanders(RTT_fwd, left_RTT_len, ref_len, alt_len,
                      transcript_strand, PAM_strand, frame_of_RTT_start,
                      RTT_genomic_positions=None, frame_map=None, boundaries=None,
                      window_nt=5, max_muts=2, max_candidates=20,
                      splice_buffer=3):
    """
    Finds synonymous ("silent") mutations that can be carried in the RTT of a
    pegRNA alongside the intended edit.

    Returns a list of dicts, each describing one bystander option:
    {'RTT': <new RTT, PAM-strand orientation>, 'n_muts': int,
     'positions': [offsets into RTT_fwd], 'dist_to_edit': int, 'aa_check': True}

    Options are sorted fewest-mutations-first, then by proximity to the edit.
    Every option is translated and checked to be synonymous against the supplied
    reading frame; all changes lie strictly within window_nt of the edit; the
    intended edit is never altered (for a pure deletion, the bases flanking the
    junction are protected instead); and positions outside the CDS or close to a
    splice site are left alone when frame_map/boundaries are provided.

    Reading frame and strand
    -------------------------
    RTT_fwd is in the PAM-strand orientation used inside
    prime.pegRNA_generator() (i.e. before reverse complementation), where:

        RTT_fwd[0:3]                  the 3 nt immediately 3' of the nick
        RTT_fwd[3:3+len(PAM)]         the PAM
        RTT_fwd[left_RTT_len:]        the intended edit

    When the transcript is on the other strand, codons are read on the reverse
    complement and the frame anchor is translated with reverse_frame_anchor().

    Indels
    -------
    Downstream of an insertion or deletion the reading frame is shifted by
    (alt_len - ref_len). Where that shift is not a multiple of three the
    downstream reading frame is destroyed and "synonymous" has no meaning, so
    codons downstream of the edit are excluded. Codons upstream of the edit are
    always safe, which is why indels are supported at all.

    Parameters
    -----------
    RTT_fwd
        *type = str*

        The RTT in PAM-strand orientation (before reverse complement), including
        the intended edit.

    left_RTT_len
        *type = int*

        Length of the left homology arm, i.e. the offset of the edit within
        RTT_fwd. Equals prime's Distance_to_nick, and excludes the edit itself.

    ref_len
        *type = int*

        Length of the reference allele that the edit replaces.

    alt_len
        *type = int*

        Length of the alternate allele carried in the RTT.

    transcript_strand
        *type = str*

        '+' or '-'; strand the transcript is on.

    PAM_strand
        *type = str*

        '+' or '-'; strand the pegRNA's PAM is on.

    frame_of_RTT_start
        *type = int*

        Codon phase (0, 1 or 2) of RTT_fwd[0], read in the transcript's
        orientation.

    RTT_genomic_positions
        *type = list or None*

        Genomic coordinate of each base of RTT_fwd, used for CDS and splice-site
        checks. Default = None (checks skipped).

    frame_map
        *type = dict or None*

        Frame map from cds_frame_map(). Required for CDS checks.

    boundaries
        *type = set or None*

        Exon boundaries from cds_boundaries(). Required for splice-site checks.

    window_nt
        *type = int*

        How far from the edit, in nt, silent mutations may be placed. Default = 5.
        This is enforced strictly: codons are pulled into the search whenever
        they overlap the window, but an option is discarded unless every base it
        changes lies inside the window itself.

    max_muts
        *type = int*

        Maximum number of silent base changes carried by a single option.
        Default = 2. Each increment multiplies the number of options roughly
        four-fold, and more changes mean more synthesis risk and more chance of
        perturbing features that are not modelled here (splice enhancers, RNA
        secondary structure).

    max_candidates
        *type = int or None*

        Maximum number of options to return, after sorting fewest-mutations
        first. Default = 20. Keeps the candidate set small enough to be scored
        downstream; set to None for no cap.

    splice_buffer
        *type = int*

        Minimum distance to keep from an exon boundary. Default = 3.
    """
    L = len(RTT_fwd)
    same_strand = (transcript_strand == PAM_strand)

    #--- put the sequence and the frame anchor in the transcript's orientation
    if same_strand:
        work_seq = RTT_fwd
        frame_offset = frame_of_RTT_start
    else:
        work_seq = str(Bio.Seq.Seq(RTT_fwd).reverse_complement())
        frame_offset = reverse_frame_anchor(frame_of_RTT_start, L)

    def to_work(offset, length=1):
        """Offset in RTT_fwd -> offset in work_seq."""
        if same_strand:
            return offset
        return L - offset - length

    def from_work(offset, length=1):
        """Offset in work_seq -> offset in RTT_fwd."""
        if same_strand:
            return offset
        return L - offset - length

    #--- the intended edit is off limits
    edit_start_fwd = left_RTT_len
    edit_end_fwd = min(left_RTT_len + alt_len, L)
    protected = set()
    for i in range(edit_start_fwd, edit_end_fwd):
        protected.add(to_work(i))

    #a pure deletion carries no alternate bases (verified against both of PEGG's
    #parsers: primedesign_formatter and mut_formatter both yield alt_seq = '' for
    #a DEL, so alt_len == 0) and would otherwise protect nothing at all. Protect
    #the bases flanking the junction instead, so that a bystander cannot be
    #placed on top of the deletion seam.
    if alt_len == 0:
        for i in (edit_start_fwd - 1, edit_start_fwd):
            if 0 <= i < L:
                protected.add(to_work(i))

    #Single reference interval for "where the edit is", used by both the window
    #and the distance metric below. Deriving those separately from
    #left_RTT_len/alt_len let their conventions drift apart -- a pure deletion
    #got a window one nt narrower on one side than a substitution did.
    if alt_len == 0:
        anchor_lo = edit_start_fwd - 1      #the two bases either side of the seam
        anchor_hi = edit_start_fwd + 1
    else:
        anchor_lo = edit_start_fwd
        anchor_hi = edit_end_fwd

    #--- frameshift: codons downstream of the edit have no defined frame
    frameshift = ((alt_len - ref_len) % 3 != 0)

    #in work_seq coordinates, is the edit before or after the rest of the RTT?
    if same_strand:
        downstream_lo = edit_end_fwd          #everything at/after this is downstream
        downstream_hi = L
    else:
        #work_seq runs along the transcript, i.e. backwards along RTT_fwd, so the
        #transcriptional downstream of the edit is the HIGH end of work_seq
        downstream_lo = to_work(edit_start_fwd) + 1
        downstream_hi = L

    #--- codon boundaries in work_seq
    frame_start_offset = (-frame_offset) % 3

    #--- the search window, in work_seq coordinates
    if same_strand:
        w_lo = max(0, anchor_lo - window_nt)
        w_hi = min(L, anchor_hi + window_nt)
    else:
        w_lo = max(0, to_work(anchor_hi - 1) - window_nt)
        w_hi = min(L, to_work(anchor_lo) + 1 + window_nt)

    #--- which codons may be varied
    first_codon = w_lo - ((w_lo - frame_start_offset) % 3)
    codon_starts = []
    c = first_codon
    while c < w_hi:
        if c >= 0 and c + 3 <= L:
            codon_starts.append(c)
        c += 3

    def codon_allowed(cs):
        """Whether the codon starting at cs in work_seq may be varied."""
        for i in range(cs, cs + 3):
            #never touch the edit
            if i in protected:
                return False
            #a frameshifting indel destroys the frame downstream of the edit
            if frameshift and downstream_lo <= i < downstream_hi:
                return False
            #CDS membership and splice-site distance
            if RTT_genomic_positions is not None and frame_map is not None:
                g = RTT_genomic_positions[from_work(i)]
                if not position_is_safe(g, frame_map, boundaries or set(),
                                        splice_buffer):
                    return False
        return True

    codon_starts = [c for c in codon_starts if codon_allowed(c)]

    if len(codon_starts) == 0:
        return []

    #span covered by the varied codons; used for the translation check below
    trans_lo = codon_starts[0]
    trans_hi = codon_starts[-1] + 3

    #--- enumerate synonymous replacements
    per_codon_options = []
    for cs in codon_starts:
        ref_codon = work_seq[cs:cs + 3].upper()
        per_codon_options.append([ref_codon] + synonymous_codons(ref_codon))

    options = []
    seen = set()

    for combo in product(*per_codon_options):
        new_work = work_seq
        for cs, codon in zip(codon_starts, combo):
            new_work = new_work[:cs] + codon + new_work[cs + 3:]

        changed = [i for i in range(L) if work_seq[i].upper() != new_work[i].upper()]

        if len(changed) == 0 or len(changed) > max_muts:
            continue
        if any(i in protected for i in changed):
            continue

        #keep the changes strictly inside the requested window. Codons are pulled
        #in whenever they overlap the window, so without this an option that only
        #touches the overhanging tail of such a codon would be returned even
        #though it sits outside window_nt of the edit.
        if any(i < w_lo or i >= w_hi for i in changed):
            continue

        #Verify synonymity rather than trusting the codon enumeration. Note what
        #this can and cannot catch: it re-translates in the frame the function
        #was given, so it will catch a mistake in the substitution or masking
        #logic, but NOT a wrong frame anchor -- a wrong anchor is self-consistent,
        #producing changes that are silent in the wrong frame and non-synonymous
        #in the real one. Only the caller can get the anchor right; see
        #reverse_frame_anchor() and the tests in tests/test_bystander.py.
        if str(Bio.Seq.Seq(work_seq[trans_lo:trans_hi]).translate()) != \
           str(Bio.Seq.Seq(new_work[trans_lo:trans_hi]).translate()):
            continue

        new_RTT = new_work if same_strand else str(
            Bio.Seq.Seq(new_work).reverse_complement())

        if new_RTT in seen:
            continue
        seen.add(new_RTT)

        positions_fwd = sorted(from_work(i) for i in changed)

        #distance from the closest silent change to the edit; mismatches close to
        #the edit are the ones that matter for MMR evasion
        dist_to_edit = min(
            min(abs(p - anchor_lo), abs(p - (anchor_hi - 1)))
            for p in positions_fwd)

        options.append({
            'RTT': new_RTT,
            'n_muts': len(changed),
            'positions': positions_fwd,
            'dist_to_edit': dist_to_edit,
            'aa_check': True,
        })

    #fewest changes first, then closest to the edit. Sorting on the sequence
    #alone would order by base identity (A < C < G < T), which biases any
    #downstream truncation towards particular bases rather than sampling evenly;
    #the sequence is kept only as a final tie-break for determinism.
    options.sort(key=lambda x: (x['n_muts'], x['dist_to_edit'], x['RTT']))

    if max_candidates is not None:
        options = options[:max_candidates]

    return options


#--- reading frame sources ---------

#Reasons a variant can end up with no silent bystanders. Recorded in the output
#so that an empty result is always explainable.
NO_FRAME_NO_ANNOTATION = 'no_frame_annotation_supplied'
NO_FRAME_OUTSIDE_CDS = 'edit_outside_supplied_CDS'
NO_FRAME_NOT_IN_FRAME = 'input_sequence_not_in_frame'


def resolve_frame_source(input_format, transcript_strand=None, start_end_cds=None,
                         ORF_start=None):
    """
    Works out where the reading frame is going to come from, and refuses the
    combinations that cannot be checked. Returns
    (mode, frame_map, boundaries) where mode is one of:

        'cds'    -- genomic coordinates plus a CDS annotation. The reading frame
                    is looked up per position, and CDS membership and splice-site
                    distance are both enforced.

        'orf'    -- no genomic coordinates. The reading frame comes from
                    ORF_start, and the caller is asserting that the whole input
                    sequence is in-frame coding sequence. CDS membership and
                    splice sites CANNOT be checked in this mode.

    Raises ValueError when silent bystanders have been asked for but the inputs
    cannot support them, rather than falling back to a guess: a wrong reading
    frame produces changes that look silent but are not, and that error would
    only surface after the library had been synthesised and sequenced.

    Parameters
    -----------
    input_format
        *type = str*

        One of 'cBioPortal', 'WT_ALT', 'PrimeDesign'.

    transcript_strand
        *type = str or None*

        '+' or '-'; the strand the transcript is on. Required for 'cBioPortal'.

    start_end_cds
        *type = list or None*

        CDS blocks for the transcript, 1-based inclusive, ordered in the +
        strand orientation. Required for 'cBioPortal'.

    ORF_start
        *type = int or None*

        0, 1 or 2: the offset within the input sequence at which the reading
        frame begins. Required for 'WT_ALT' and 'PrimeDesign'.
    """
    if input_format == 'cBioPortal':
        if start_end_cds is None or transcript_strand is None:
            raise ValueError(
                "silent_bystander=True with input_format='cBioPortal' requires "
                "both start_end_cds and transcript_strand, so that the reading "
                "frame can be looked up per genomic position. Supply them, or "
                "set silent_bystander=False.")

        frame_map = cds_frame_map(start_end_cds, transcript_strand)
        boundaries = cds_boundaries(start_end_cds)
        return 'cds', frame_map, boundaries

    #WT_ALT / PrimeDesign: no genomic coordinates are available
    if ORF_start is None:
        raise ValueError(
            "silent_bystander=True with input_format='%s' requires ORF_start "
            "(0, 1 or 2).\n"
            "These formats carry no genomic coordinates, so the reading frame "
            "cannot be looked up and must be declared: the input sequence has "
            "to be IN FRAME, i.e. entirely coding sequence with the frame "
            "beginning at offset ORF_start.\n"
            "Note that CDS membership and splice sites cannot be checked in "
            "this mode. If any part of the input sequence is intronic or "
            "untranslated, use input_format='cBioPortal' with start_end_cds "
            "instead, or set silent_bystander=False." % input_format)

    if ORF_start not in [0, 1, 2]:
        raise ValueError("ORF_start must be 0, 1 or 2 (got %r)" % (ORF_start,))

    return 'orf', None, None


def frame_of_offset(offset, mode, genomic_positions=None, frame_map=None,
                    ORF_start=None):
    """
    Returns the codon phase (0, 1 or 2) at a given offset of the design
    sequence, or None if it cannot be determined.

    Parameters
    -----------
    offset
        *type = int*

        Offset within the design sequence.

    mode
        *type = str*

        'cds' or 'orf'; see resolve_frame_source().

    genomic_positions
        *type = list or None*

        Genomic coordinate of every base of the design sequence. Required for
        'cds' mode. See prime.mutation.genomic_positions().

    frame_map
        *type = dict or None*

        Frame map from cds_frame_map(). Required for 'cds' mode.

    ORF_start
        *type = int or None*

        Frame offset within the design sequence. Required for 'orf' mode.
    """
    if mode == 'cds':
        if genomic_positions is None or frame_map is None:
            return None
        entry = frame_at(genomic_positions[offset], frame_map)
        if entry is None:
            return None
        return entry[1]

    #'orf': the frame is declared relative to the start of the input sequence
    return (offset - ORF_start) % 3
