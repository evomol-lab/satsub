"""Alignment input/output utilities.

Handles reading alignments from common formats (FASTA, Clustal, Phylip,
Nexus, Stockholm), detecting whether a file holds nucleotide or amino-acid
sequences, normalizing RNA to DNA for internal processing, and translating
in-frame nucleotide alignments to amino acids.
"""

from __future__ import annotations

import io
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

SeqKind = Literal["nucleotide", "amino_acid", "unknown"]

# Formats attempted, in order, when the format is not explicitly given.
_CANDIDATE_FORMATS = ["fasta", "clustal", "phylip-relaxed", "phylip", "nexus", "stockholm"]

# Characters that legitimately appear in nucleotide alignments (DNA or RNA),
# including IUPAC ambiguity codes, gap and missing-data symbols.
_NUCLEOTIDE_ALPHABET = set("ACGTUN RYSWKMBDHV-.?".replace(" ", ""))

# The 20 standard amino-acid one-letter codes plus common ambiguity codes.
_AMINO_ACID_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYBZJXU*O-.?")


class AlignmentReadError(ValueError):
    """Raised when an alignment file cannot be parsed in any known format."""


@dataclass
class LoadedAlignment:
    alignment: MultipleSeqAlignment
    kind: SeqKind
    format: str
    source_name: str
    had_rna: bool = False
    warnings: list[str] = field(default_factory=list)

    @property
    def n_sequences(self) -> int:
        return len(self.alignment)

    @property
    def length(self) -> int:
        return self.alignment.get_alignment_length()

    @property
    def names(self) -> list[str]:
        return [rec.id for rec in self.alignment]


def _read_text(file_like_or_path) -> tuple[str, str]:
    """Return (text, source_name) from a path, str/bytes path, or a
    file-like/UploadedFile object (as produced by Streamlit's uploader)."""
    if hasattr(file_like_or_path, "read"):
        raw = file_like_or_path.read()
        name = getattr(file_like_or_path, "name", "uploaded_file")
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8", errors="replace")
        return raw, name
    path = Path(file_like_or_path)
    return path.read_text(errors="replace"), path.name


def sniff_and_parse(text: str) -> tuple[MultipleSeqAlignment, str]:
    """Try each supported alignment format until one parses successfully.

    Returns the parsed alignment and the format name that worked.
    """
    last_error: Exception | None = None
    for fmt in _CANDIDATE_FORMATS:
        try:
            aln = AlignIO.read(io.StringIO(text), fmt)
            if len(aln) >= 1 and aln.get_alignment_length() >= 1:
                return aln, fmt
        except Exception as exc:  # noqa: BLE001 - deliberately broad, we try many formats
            last_error = exc
            continue
    raise AlignmentReadError(
        "Could not parse alignment in any supported format "
        f"(tried: {', '.join(_CANDIDATE_FORMATS)}). Last error: {last_error}"
    )


def detect_kind(alignment: MultipleSeqAlignment, sample_size: int = 20) -> SeqKind:
    """Classify an alignment as nucleotide or amino acid by character content.

    Looks at up to `sample_size` sequences and decides based on what fraction
    of non-gap characters fall within the nucleotide alphabet. Protein
    sequences that happen to be rich in A/C/G/T/N letters are still correctly
    classified because the amino-acid alphabet is a superset check performed
    second, and pure-ACGTUN content is overwhelmingly likely to be nucleic
    acid in practice.
    """
    letters = []
    for rec in list(alignment)[:sample_size]:
        letters.extend(str(rec.seq).upper())
    non_gap = [c for c in letters if c not in "-.?"]
    if not non_gap:
        return "unknown"
    nuc_count = sum(1 for c in non_gap if c in _NUCLEOTIDE_ALPHABET)
    frac_nuc = nuc_count / len(non_gap)
    # Nucleotide alphabet (ACGTUN + ambiguity codes) is a subset of letters
    # also valid as amino acid codes, so we require a high fraction AND that
    # the "core" unambiguous bases (ACGTU) dominate, not just N/ambiguity codes.
    core = set("ACGTU")
    frac_core = sum(1 for c in non_gap if c in core) / len(non_gap)
    if frac_nuc >= 0.95 and frac_core >= 0.75:
        return "nucleotide"
    return "amino_acid"


def has_rna(alignment: MultipleSeqAlignment) -> bool:
    for rec in alignment:
        if "U" in str(rec.seq).upper():
            return True
    return False


def normalize_nucleotide(alignment: MultipleSeqAlignment) -> MultipleSeqAlignment:
    """Uppercase all sequences and convert RNA (U) to DNA (T) in place semantics
    (returns a new alignment object; does not mutate the input)."""
    records = []
    for rec in alignment:
        seq_str = str(rec.seq).upper().replace("U", "T")
        new_rec = SeqRecord(Seq(seq_str), id=rec.id, name=rec.name, description=rec.description)
        records.append(new_rec)
    return MultipleSeqAlignment(records)


def read_alignment(file_like_or_path) -> LoadedAlignment:
    """Read an alignment from a path or a file-like object and classify it.

    Nucleotide alignments are normalized (uppercase, U->T) before being
    returned; the original had_rna flag records whether the source used the
    RNA alphabet so the GUI/reports can mention it.
    """
    text, source_name = _read_text(file_like_or_path)
    alignment, fmt = sniff_and_parse(text)

    warnings: list[str] = []
    lengths = {len(rec.seq) for rec in alignment}
    if len(lengths) > 1:
        warnings.append(
            "Sequences have unequal lengths after parsing; this file may not "
            "be a proper alignment (no gap-padding to a common length)."
        )

    kind = detect_kind(alignment)
    rna = has_rna(alignment) if kind == "nucleotide" else False
    if kind == "nucleotide":
        alignment = normalize_nucleotide(alignment)
    else:
        # Still uppercase amino acid sequences for consistent comparisons.
        records = [
            SeqRecord(Seq(str(rec.seq).upper()), id=rec.id, name=rec.name, description=rec.description)
            for rec in alignment
        ]
        alignment = MultipleSeqAlignment(records)

    # Deduplicate/repair empty or duplicate sequence names, which are common
    # in real-world FASTA files and would otherwise break pairwise reporting.
    seen: dict[str, int] = {}
    for rec in alignment:
        name = rec.id or "seq"
        if name in seen:
            seen[name] += 1
            warnings.append(f"Duplicate sequence name '{name}' encountered; disambiguating.")
            rec.id = f"{name}_{seen[name]}"
        else:
            seen[name] = 0

    return LoadedAlignment(
        alignment=alignment,
        kind=kind,
        format=fmt,
        source_name=source_name,
        had_rna=rna,
        warnings=warnings,
    )


# Standard genetic code tables commonly needed for barcoding/coding-sequence
# work, keyed by the NCBI translation table id used by Biopython.
GENETIC_CODE_CHOICES: dict[int, str] = {
    1: "Standard",
    2: "Vertebrate Mitochondrial",
    3: "Yeast Mitochondrial",
    4: "Mold/Protozoan/Coelenterate Mitochondrial",
    5: "Invertebrate Mitochondrial",
    6: "Ciliate/Dasycladacean/Hexamita Nuclear",
    9: "Echinoderm/Flatworm Mitochondrial",
    10: "Euplotid Nuclear",
    11: "Bacterial/Archaeal/Plant Plastid",
    12: "Alternative Yeast Nuclear",
    13: "Ascidian Mitochondrial",
    14: "Alternative Flatworm Mitochondrial",
    16: "Chlorophycean Mitochondrial",
    21: "Trematode Mitochondrial",
    22: "Scenedesmus obliquus Mitochondrial",
    23: "Thraustochytrium Mitochondrial",
    24: "Pterobranchia Mitochondrial",
    25: "Candidate Division SR1 and Gracilibacteria",
}


def translate_alignment(
    alignment: MultipleSeqAlignment,
    frame_start: int = 1,
    table: int = 1,
) -> tuple[MultipleSeqAlignment, list[str]]:
    """Translate an in-frame nucleotide alignment to amino acids.

    Columns are consumed in triplets starting at `frame_start` (1-based: 1,
    2 or 3). A trailing partial codon (fewer than 3 columns remaining) is
    dropped. A codon containing any gap character is translated to '-' if it
    is entirely gaps, or 'X' if it is a mix of gap and non-gap positions
    (ambiguous / frameshift-affected). Stop codons are translated to '*'.

    Returns the translated alignment and a list of human-readable warnings
    (e.g. sequences containing internal stop codons, which often signals the
    wrong genetic code table or reading frame was chosen).
    """
    warnings: list[str] = []
    start = frame_start - 1
    usable_len = ((alignment.get_alignment_length() - start) // 3) * 3
    if usable_len <= 0:
        raise ValueError("Alignment too short for the selected reading frame.")

    records = []
    for rec in alignment:
        nt = str(rec.seq).upper()[start : start + usable_len]
        aa_chars = []
        has_internal_stop = False
        n_codons = usable_len // 3
        for i in range(n_codons):
            codon = nt[i * 3 : i * 3 + 3]
            if codon == "---":
                aa_chars.append("-")
                continue
            if "-" in codon or any(c not in "ACGT" for c in codon):
                aa_chars.append("X")
                continue
            aa = str(Seq(codon).translate(table=table))
            if aa == "*" and i < n_codons - 1:
                has_internal_stop = True
            aa_chars.append(aa)
        if has_internal_stop:
            warnings.append(
                f"Sequence '{rec.id}' has an internal stop codon under table {table} "
                f"starting at frame {frame_start}; double-check the genetic code and reading frame."
            )
        records.append(SeqRecord(Seq("".join(aa_chars)), id=rec.id, name=rec.name, description=""))

    return MultipleSeqAlignment(records), warnings


def example_alignment_path() -> Path:
    """Path to the bundled Vertebrate COI example alignment."""
    return Path(__file__).parent / "example_data" / "VertCOI.fas"


def alignment_from_pairs(pairs) -> MultipleSeqAlignment:
    """Rebuild an alignment from a sequence of (name, sequence_string) pairs.

    Used together with `alignment_to_pairs` to pass alignments through
    Streamlit's hash-based caching using a plain, hashable representation.
    """
    records = [SeqRecord(Seq(seq), id=name, name=name, description="") for name, seq in pairs]
    return MultipleSeqAlignment(records)


def alignment_to_pairs(alignment: MultipleSeqAlignment) -> tuple[tuple[str, str], ...]:
    """Inverse of `alignment_from_pairs`: a hashable snapshot of an alignment."""
    return tuple((rec.id, str(rec.seq)) for rec in alignment)
