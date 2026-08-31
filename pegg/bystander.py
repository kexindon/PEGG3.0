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
equally, because they do not all carry genomic coordinates::

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

    Returns a list of dicts, each describing one bystander option, with keys
    'RTT' (the new RTT in PAM-strand orientation), 'n_muts', 'positions'
    (offsets into RTT_fwd), 'dist_to_edit' and 'aa_check'.

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


#--- building CDS annotations ---------

def cds_from_gtf(filepath, genes=None, transcript_ids=None, longest_only=False):
    """
    Extracts CDS blocks per gene from a GTF/GFF annotation file, in the format
    prime.run() expects. Returns a dict of

        {gene_name: {'strand': '+'/'-', 'cds': [[start, end], ...],
                     'transcript_id': str, 'n_codons': int, 'valid': bool}}

    Blocks are 1-based inclusive and ordered along the + strand, matching
    library.neutral_substitutions() and prime.run(start_end_cds=...).

    A gene typically has several transcripts, and the reading frame belongs to
    one of them, so exactly one has to be chosen. By default the transcript with
    the longest CDS is used; pass transcript_ids to choose explicitly.

    Note the chromosome names in the annotation must match the keys of the
    chrom_dict used for design (PEGG uses 1..22, 'X', 'Y'); see
    cds_annotation_report() for a check of this and of frame validity.

    Parameters
    -----------
    filepath
        *type = str*

        Path to a GTF or GFF3 file, optionally gzipped. Ensembl, GENCODE and
        NCBI RefSeq files all work.

    genes
        *type = list or None*

        Gene names to extract. Default = None (every gene in the file, which is
        slow and memory-hungry for a whole genome).

    transcript_ids
        *type = dict or None*

        {gene_name: transcript_id} to pin specific transcripts instead of taking
        the longest CDS. Default = None.

    longest_only
        *type = bool*

        Whether to reduce each gene to its longest-CDS transcript. Default =
        False, which returns every transcript keyed as "gene|transcript_id" so
        that the choice stays explicit -- the reading frame belongs to one
        transcript, and picking the longest is a guess that can disagree with
        the canonical transcript a variant was annotated against.

        Pass transcript_ids={gene: id} to select transcripts by name; set this
        to True only if an arbitrary-but-reproducible choice is acceptable.
    """
    import gzip
    import re as _re

    opener = gzip.open if str(filepath).endswith('.gz') else open

    wanted = set(genes) if genes is not None else None

    #transcript -> (gene, strand, chrom, [[start, end], ...])
    transcripts = {}

    with opener(filepath, 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                continue

            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9 or parts[2] != 'CDS':
                continue

            chrom, start, end, strand, attrs = (parts[0], int(parts[3]),
                                                int(parts[4]), parts[6], parts[8])

            #GTF uses key "value"; GFF3 uses key=value
            gene = (_re.search(r'gene_name[ =]"?([^";]+)"?', attrs)
                    or _re.search(r'\bgene[ =]"?([^";]+)"?', attrs))
            tx = (_re.search(r'transcript_id[ =]"?([^";]+)"?', attrs)
                  or _re.search(r'\bParent[ =]"?(?:rna[-:])?([^";]+)"?', attrs))

            if gene is None or tx is None:
                continue

            gene = gene.group(1)
            if wanted is not None and gene not in wanted:
                continue

            tx = tx.group(1)
            entry = transcripts.setdefault(tx, {'gene': gene, 'strand': strand,
                                                'chrom': chrom, 'cds': []})
            entry['cds'].append([start, end])

    #pick one transcript per gene
    by_gene = {}
    for tx, entry in transcripts.items():
        by_gene.setdefault(entry['gene'], []).append((tx, entry))

    result = {}
    for gene, candidates in by_gene.items():
        if transcript_ids is not None and gene in transcript_ids:
            chosen = [(tx, e) for tx, e in candidates
                      if tx == transcript_ids[gene] or
                      tx.split('.')[0] == str(transcript_ids[gene]).split('.')[0]]
            if len(chosen) == 0:
                continue
            candidates = chosen
        elif longest_only:
            candidates = [max(candidates,
                              key=lambda x: sum(e - s + 1 for s, e in x[1]['cds']))]

        for tx, entry in candidates:
            blocks = sorted(entry['cds'])
            total = sum(e - s + 1 for s, e in blocks)
            key = gene if (longest_only and transcript_ids is None) or len(candidates) == 1 \
                else '%s|%s' % (gene, tx)
            result[key] = {'strand': entry['strand'],
                           'cds': blocks,
                           'chrom': entry['chrom'],
                           'transcript_id': tx,
                           'n_codons': total // 3,
                           'valid': (total % 3 == 0) and total > 0}

    return result


def cds_annotation_report(gene_cds, mutations=None, gene_column='Hugo_Symbol'):
    """
    Checks a CDS annotation dict before it is used for design, and returns a
    dataframe summarising it. Catches the mistakes that would otherwise only
    surface as silently wrong bystanders:

      - a CDS whose length is not a multiple of three (wrong or partial transcript)
      - overlapping CDS blocks
      - a gene in the variant table with no annotation (gets no bystanders)
      - variants that fall outside the CDS of their own gene (also get none)

    Parameters
    -----------
    gene_cds
        *type = dict*

        {gene: {'strand': ..., 'cds': [[start, end], ...]}}, as produced by
        cds_from_gtf() or written by hand.

    mutations
        *type = pd.DataFrame or None*

        The input variant table, to cross-check gene coverage and whether each
        variant lies inside its gene's CDS. Default = None (skip that check).

    gene_column
        *type = str*

        Column of the variant table holding the gene name.
        Default = 'Hugo_Symbol'.
    """
    import pandas as pd

    rows = []
    for gene, ann in gene_cds.items():
        blocks = sorted(ann['cds'])
        total = sum(e - s + 1 for s, e in blocks)

        overlaps = any(blocks[i][1] >= blocks[i + 1][0]
                       for i in range(len(blocks) - 1))

        problems = []
        if total == 0:
            problems.append('empty CDS')
        if total % 3 != 0:
            problems.append('not a multiple of 3')
        if overlaps:
            problems.append('overlapping blocks')
        if ann.get('strand') not in ('+', '-'):
            problems.append('strand must be + or -')

        row = {'gene': gene,
               'strand': ann.get('strand'),
               'n_blocks': len(blocks),
               'cds_bp': total,
               'n_codons': total / 3,
               'transcript_id': ann.get('transcript_id'),
               'n_variants': 0,
               'variants_in_cds': 0,
               'problems': '; '.join(problems) if problems else ''}

        if mutations is not None:
            sub = mutations[mutations[gene_column] == gene]
            row['n_variants'] = len(sub)
            if len(sub) and not problems:
                fmap = cds_frame_map(blocks, ann['strand'])
                row['variants_in_cds'] = int(sum(
                    1 for p in sub['Start_Position'] if int(p) in fmap))
                if row['variants_in_cds'] < len(sub):
                    missing = len(sub) - row['variants_in_cds']
                    row['problems'] = ('%d/%d variants outside the CDS'
                                       % (missing, len(sub)))

        rows.append(row)

    report = pd.DataFrame(rows)

    #genes present in the variant table but not annotated at all
    if mutations is not None:
        missing = sorted(set(mutations[gene_column]) - set(gene_cds))
        for gene in missing:
            report = pd.concat([report, pd.DataFrame([{
                'gene': gene, 'strand': None, 'n_blocks': 0, 'cds_bp': 0,
                'n_codons': 0, 'transcript_id': None,
                'n_variants': int((mutations[gene_column] == gene).sum()),
                'variants_in_cds': 0,
                'problems': 'no CDS annotation - will get no bystanders'}])],
                ignore_index=True)

    return report.sort_values('gene').reset_index(drop=True)


def cds_from_annotation_db(db, gene_or_tx):
    """
    Extracts CDS blocks and strand from a gffutils annotation database, e.g. the
    one H2M's anno_loader() returns.

    Returns a dict with keys 'strand', 'cds', 'chrom', 'transcript_id',
    'n_codons' and 'valid' -- the same shape as cds_from_gtf(), so it drops
    straight into prime.run().

    The CDS extraction follows H2M's own GetTx(): the transcript's CDS children
    are pulled in + strand order and used as the coding blocks. Reimplemented
    here rather than imported so that pegg keeps no dependency on H2M; gffutils
    is the only requirement, and only when this function is called.

    Parameters
    -----------
    db
        *type = gffutils.FeatureDB*

        Annotation database, e.g. from gffutils.FeatureDB(path) or
        h2m.anno_loader(path). Must be the same genome build as the chrom_dict
        used for design.

    gene_or_tx
        *type = str*

        A transcript id (e.g. 'ENST00000269305'), or a gene name / gene id, in
        which case the transcript with the longest CDS is used.

    """
    import gffutils

    def _blocks(tx_id):
        #H2M's GetTx pulls CDS children ordered along the + strand
        cds = list(db.children(tx_id, order_by='+end', featuretype=['CDS']))
        return [[int(i.start), int(i.end)] for i in cds]

    #resolve gene name / gene id to a transcript
    tx_id = gene_or_tx
    feature = None
    is_tx = False

    #Transcript ids may or may not carry a version suffix (ENST00000269305.4 vs
    #ENST00000269305), and the annotation may use either form, so try both.
    bare = str(gene_or_tx).split('.')[0]
    for candidate in ([gene_or_tx] if bare == gene_or_tx else [gene_or_tx, bare]):
        try:
            feature = db[candidate]
        except gffutils.FeatureNotFoundError:
            continue
        if feature.featuretype in ('transcript', 'mRNA'):
            tx_id, is_tx = candidate, True
            break

    #still nothing, but the id looks like a versioned transcript: scan for a
    #transcript whose id matches once the version is stripped
    if not is_tx and bare != gene_or_tx:
        for t in db.features_of_type(('transcript', 'mRNA')):
            if str(t.id).split('.')[0] == bare:
                tx_id, feature, is_tx = t.id, t, True
                break

    if not is_tx:
        #A gene has many transcripts and the reading frame belongs to exactly one
        #of them, so picking one automatically would be a guess. Guessing wrong
        #does not fail loudly: it yields bystanders that are silent in the wrong
        #frame, and protein-level annotation (e.g. HGVSp_Short) that no longer
        #lines up with the transcript actually designed against. Make the caller
        #choose.
        candidates = []

        if feature is not None and feature.featuretype == 'gene':
            candidates = [t.id for t in db.children(gene_or_tx,
                                                    featuretype=['transcript', 'mRNA'])]
        else:
            for g in db.features_of_type('gene'):
                names = g.attributes.get('gene_name', []) + g.attributes.get('gene', [])
                if gene_or_tx in names:
                    candidates = [t.id for t in db.children(g.id,
                                                            featuretype=['transcript', 'mRNA'])]
                    break

        if len(candidates) == 0:
            raise ValueError("no transcript found for %r in this annotation "
                             "database" % (gene_or_tx,))

        raise ValueError(
            "%r is a gene, not a transcript, and the reading frame belongs to a "
            "single transcript -- pass a transcript id instead of letting it be "
            "guessed.\n"
            "This annotation has %d transcript(s) for it, e.g. %s.\n"
            "Use the canonical transcript for your variant annotation: H2M's "
            "get_tx_batch(df, species, ver) attaches canonical ids to a variant "
            "table, or pass transcript_ids={%r: '<id>'}."
            % (gene_or_tx, len(candidates), ', '.join(sorted(candidates)[:3]),
               gene_or_tx))

    blocks = sorted(_blocks(tx_id))
    if len(blocks) == 0:
        raise ValueError("transcript %r has no CDS features (non-coding?)" % (tx_id,))

    total = sum(e - s + 1 for s, e in blocks)
    tx = db[tx_id]

    return {'strand': tx.strand,
            'cds': blocks,
            'chrom': tx.chrom,
            'transcript_id': tx_id,
            'n_codons': total // 3,
            'valid': (total % 3 == 0) and total > 0}


def cds_from_tx_ids(db, tx_ids):
    """
    Looks up CDS blocks and strand for a set of transcript ids at once. Returns
    {tx_id: {'strand': ..., 'cds': [[start, end], ...], ...}}, with one entry per
    transcript that could be resolved; transcripts that are missing from the
    annotation or have no CDS are skipped rather than raising.

    Designed to take the transcript ids that H2M's get_tx_batch() attaches to a
    variant table, so that a library spanning many genes can be designed without
    naming transcripts one at a time.

    Parameters
    -----------
    db
        *type = gffutils.FeatureDB*

        Annotation database, e.g. gffutils.FeatureDB(path) or the object
        h2m.anno_loader(path) returns.

    tx_ids
        *type = iterable of str*

        Transcript ids, e.g. list(df['tx_id_h'].unique()).
    """
    result = {}

    for tx in tx_ids:
        if tx is None or (isinstance(tx, float) and tx != tx):   #NaN
            continue
        if tx in result:
            continue
        try:
            result[tx] = cds_from_annotation_db(db, tx)
        except (ValueError, KeyError):
            continue

    return result


def add_cds_to_variants(df, db, tx_column='tx_id_h', gene_column=None):
    """
    Attaches reading frame information to a variant table, one transcript per
    row, and returns (annotated_df, cds_lookup).

    The returned dataframe gains 'transcript_strand', 'cds_valid' and
    'cds_n_codons' columns, and cds_lookup maps each transcript id to its CDS
    blocks. Together these are what prime.run(silent_bystander=True) needs, and
    they keep each variant tied to the transcript it was annotated against --
    which matters when a library spans several genes, since the reading frame
    belongs to one transcript and cannot be shared between them.

    Typical use, starting from H2M::

        df, df_fail = h2m.get_tx_batch(df, species='h', ver=37)
        df, cds_lookup = bystander.add_cds_to_variants(df, db_h)

        for tx, sub in df.groupby('tx_id_h'):
            ann = cds_lookup[tx]
            peg = prime.run(sub, 'cBioPortal', chrom_dict=chrom_dict,
                            silent_bystander=True,
                            transcript_strand=ann['strand'],
                            start_end_cds=ann['cds'], ...)

    Parameters
    -----------
    df
        *type = pd.DataFrame*

        Variant table carrying a transcript id column.

    db
        *type = gffutils.FeatureDB*

        Annotation database; must be the same genome build as the variants.

    tx_column
        *type = str*

        Column holding the transcript id. Default = 'tx_id_h' (H2M's human
        output); use 'tx_id_m' for mouse.

    gene_column
        *type = str or None*

        Gene name column, used only to make the summary printout more readable.
        Default = None.
    """
    import pandas as pd

    if tx_column not in df.keys():
        raise ValueError(
            "no %r column in the variant table. Run H2M's get_tx_batch() first, "
            "or pass tx_column= the column holding your transcript ids."
            % (tx_column,))

    cds_lookup = cds_from_tx_ids(db, df[tx_column].dropna().unique())

    out = df.copy()
    out['transcript_strand'] = [
        cds_lookup[t]['strand'] if t in cds_lookup else None
        for t in out[tx_column]]
    out['cds_valid'] = [
        bool(cds_lookup[t]['valid']) if t in cds_lookup else False
        for t in out[tx_column]]
    out['cds_n_codons'] = [
        cds_lookup[t]['n_codons'] if t in cds_lookup else None
        for t in out[tx_column]]

    #say plainly which variants will and will not get bystanders
    n_ok = int(out['cds_valid'].sum())
    print('%d/%d variants have a usable reading frame (%d transcripts)'
          % (n_ok, len(out), len(cds_lookup)))

    if n_ok < len(out):
        bad = out[~out['cds_valid']]
        cols = [c for c in (gene_column, tx_column) if c is not None]
        missing = sorted(set(bad[tx_column].dropna()) - set(cds_lookup))
        if missing:
            print('  no CDS in the annotation for: %s'
                  % ', '.join(str(m) for m in missing[:8]))
        invalid = sorted({t for t in bad[tx_column].dropna()
                          if t in cds_lookup and not cds_lookup[t]['valid']})
        if invalid:
            print('  CDS not a multiple of 3 for: %s'
                  % ', '.join(str(m) for m in invalid[:8]))
        print('  these variants will still get ordinary pegRNAs, just no bystanders')

    return out, cds_lookup


def cds_for_variants(df, db, gene_column='Hugo_Symbol', tx_column=None,
                     transcript_ids=None, species='human', genome_version=37,
                     verbose=True):
    """
    Resolves the reading frame for a variant table, working from the standard
    cBioPortal column names, and returns (annotated_df, cds_lookup).

    This is the batch entry point for cBioPortal-format input: it works from the
    gene names PEGG already uses ('Hugo_Symbol') rather than requiring H2M's
    'gene_name_h' column, while still leaving the choice of transcript explicit.

    The transcript of each gene is its canonical transcript, from the curated
    table H2M uses (see canonical_transcripts()). This matters because variant
    annotation -- HGVSp_Short, protein position -- is generally computed against
    the canonical transcript, so designing against a different one can shift the
    reading frame and make those labels disagree with the design. Override per
    gene with transcript_ids, or use ids already in the table with tx_column.

    The returned dataframe gains 'transcript_id', 'transcript_strand',
    'cds_valid' and 'cds_n_codons'; cds_lookup maps each gene to its CDS blocks.
    Genes that cannot be resolved are reported and simply get no bystanders.

    If a transcript id column is already present -- for instance because
    h2m.get_tx_batch() was run first -- pass tx_column to use those ids instead
    of choosing per gene.

    Typical use::

        df, cds_lookup = bystander.cds_for_variants(mutations, db)

        for gene, sub in df.groupby('Hugo_Symbol'):
            ann = cds_lookup.get(gene)
            if ann is not None and ann['valid']:
                peg = prime.run(sub.reset_index(drop=True), 'cBioPortal',
                                silent_bystander=True,
                                transcript_strand=ann['strand'],
                                start_end_cds=ann['cds'], ...)

    Parameters
    -----------
    df
        *type = pd.DataFrame*

        Variant table in cBioPortal format.

    db
        *type = gffutils.FeatureDB*

        Annotation database, e.g. gffutils.FeatureDB(path) or the object
        h2m.anno_loader(path) returns. Must be the same genome build as the
        variants.

    gene_column
        *type = str*

        Column holding the gene name. Default = 'Hugo_Symbol' (cBioPortal).

    tx_column
        *type = str or None*

        Column holding a transcript id, if the table already has one. H2M's
        get_tx_batch() writes canonical ids into 'tx_id_h' / 'tx_id_m'.
        Default = None, i.e. look the canonical transcript up per gene.

    species
        *type = str*

        'human' or 'mouse', for the canonical transcript table.
        Default = 'human'.

    genome_version
        *type = int*

        37 or 38, matching the annotation database and the variant coordinates.
        Default = 37.

    transcript_ids
        *type = dict or None*

        {gene: transcript_id} overriding the canonical transcript for the genes
        it names. Default = None (use the canonical transcript throughout).

    verbose
        *type = bool*

        Whether to print a summary of what resolved and what did not.
        Default = True.
    """
    import pandas as pd

    key_column = tx_column if tx_column is not None else gene_column

    if key_column not in df.keys():
        raise ValueError(
            "no %r column in the variant table. Pass gene_column= or "
            "tx_column= to point at the column holding gene names or "
            "transcript ids." % (key_column,))

    #The reading frame belongs to one transcript, so a transcript is named rather
    #than guessed. Unless the table already carries ids (tx_column), the
    #canonical transcript of each gene is used -- the same curated table H2M
    #uses -- with transcript_ids overriding it where given.
    canonical = None
    if tx_column is None:
        canonical = canonical_transcripts(species, genome_version)

    lookup = {}
    failures = {}

    for key in pd.Series(df[key_column]).dropna().unique():
        target = key
        if tx_column is None:
            if transcript_ids is not None and key in transcript_ids:
                target = transcript_ids[key]
            elif key in canonical:
                target = canonical[key]
            else:
                failures[key] = 'no canonical transcript known for this gene'
                continue

        try:
            lookup[key] = cds_from_annotation_db(db, target)
        except (ValueError, KeyError) as err:
            failures[key] = str(err).split('\n')[0]

    out = df.copy()
    out['transcript_id'] = [lookup[k]['transcript_id'] if k in lookup else None
                            for k in out[key_column]]
    out['transcript_strand'] = [lookup[k]['strand'] if k in lookup else None
                                for k in out[key_column]]
    out['cds_valid'] = [bool(lookup[k]['valid']) if k in lookup else False
                        for k in out[key_column]]
    out['cds_n_codons'] = [lookup[k]['n_codons'] if k in lookup else None
                           for k in out[key_column]]

    if verbose:
        n_ok = int(out['cds_valid'].sum())
        print('%d/%d variants have a usable reading frame (%d of %d %s resolved)'
              % (n_ok, len(out), len(lookup),
                 out[key_column].nunique(),
                 'transcripts' if tx_column else 'genes'))

        if failures:
            print('  not found in the annotation: %s'
                  % ', '.join(sorted(failures)[:8]))
        bad_frame = sorted(k for k, v in lookup.items() if not v['valid'])
        if bad_frame:
            print('  CDS not a multiple of 3: %s' % ', '.join(bad_frame[:8]))
        if failures or bad_frame:
            print('  those variants still get ordinary pegRNAs, just no bystanders')

    return out, lookup


#--- canonical transcripts ---------

#Canonical transcript per gene, taken from H2M
#(https://github.com/kexindong/h2m-public), which curates one transcript per gene
#rather than leaving the choice to a heuristic. Human entries are
#[GRCh37, GRCh38]; mouse entries are a single id.
_CANONICAL_TX_CACHE = {}


def canonical_transcripts(species='human', genome_version=37):
    """
    Returns {gene_name: transcript_id}, the canonical transcript of each gene.

    The table is the one H2M uses, bundled here so that pegg needs no dependency
    on H2M. Using the canonical transcript matters because variant annotation
    (HGVSp_Short, protein position) is generally computed against it: designing
    against a different transcript can shift the reading frame and make the
    protein-level labels disagree with what was actually designed.

    Parameters
    -----------
    species
        *type = str*

        'human' or 'mouse'. Default = 'human'.

    genome_version
        *type = int*

        37 or 38, for human. Ignored for mouse. Default = 37.
    """
    assert species in ['human', 'mouse'], "species must be 'human' or 'mouse'"
    assert genome_version in [37, 38], 'genome_version must be 37 or 38'

    key = (species, genome_version if species == 'human' else None)
    if key in _CANONICAL_TX_CACHE:
        return _CANONICAL_TX_CACHE[key]

    import json
    from importlib.resources import files

    filename = ('canonical_tx_human.json' if species == 'human'
                else 'canonical_tx_mouse.json')
    with open(files(__package__).joinpath(filename)) as handle:
        raw = json.load(handle)

    if species == 'human':
        #entries are [GRCh37 id, GRCh38 id]
        index = 0 if genome_version == 37 else 1
        table = {gene: ids[index] for gene, ids in raw.items()
                 if isinstance(ids, list) and len(ids) > index}
    else:
        table = dict(raw)

    _CANONICAL_TX_CACHE[key] = table
    return table


def canonical_transcript(gene, species='human', genome_version=37):
    """
    Returns the canonical transcript id of one gene, or None if the gene is not
    in the table.

    Parameters
    -----------
    gene
        *type = str*

        Gene name, e.g. 'TP53'.

    species
        *type = str*

        'human' or 'mouse'. Default = 'human'.

    genome_version
        *type = int*

        37 or 38, for human. Default = 37.
    """
    return canonical_transcripts(species, genome_version).get(gene)
