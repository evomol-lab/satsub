"""SatSub Streamlit GUI.

Run with:  streamlit run satsub/app.py
or, once installed:  satsub-gui
"""

from __future__ import annotations

import pandas as pd
import streamlit as st

from satsub.io import (
    GENETIC_CODE_CHOICES,
    AlignmentReadError,
    alignment_from_pairs,
    alignment_to_pairs,
    example_alignment_path,
    read_alignment,
    translate_alignment,
)
from satsub.plotting import categorical_bar, identity_heatmap, monochrome_bar, saturation_scatter
from satsub.saturation import MODEL_LABELS, POSITION_LABELS, compute_saturation_table
from satsub.stats_aa import summarize_amino_acid_alignment
from satsub.stats_nt import summarize_nucleotide_alignment

st.set_page_config(page_title="SatSub", page_icon="\U0001F9EC", layout="wide")


# --------------------------------------------------------------------------
# Cached computations (keyed on a hashable snapshot of the alignment, not
# the Biopython object itself, so results survive unrelated UI interactions)
# --------------------------------------------------------------------------


@st.cache_data(show_spinner="Computing nucleotide statistics...")
def _cached_nt_summary(pairs, frame_start):
    aln = alignment_from_pairs(pairs)
    return summarize_nucleotide_alignment(aln, frame_start=frame_start)


@st.cache_data(show_spinner="Computing amino-acid statistics...")
def _cached_aa_summary(pairs):
    aln = alignment_from_pairs(pairs)
    return summarize_amino_acid_alignment(aln)


@st.cache_data(show_spinner="Computing saturation table...")
def _cached_saturation(pairs, subset, frame_start):
    aln = alignment_from_pairs(pairs)
    return compute_saturation_table(aln, subset=subset, frame_start=frame_start)


@st.cache_data
def _cached_translate(pairs, frame_start, table):
    aln = alignment_from_pairs(pairs)
    aa_aln, warnings = translate_alignment(aln, frame_start=frame_start, table=table)
    return alignment_to_pairs(aa_aln), warnings


def _df_download_button(df: pd.DataFrame, label: str, filename: str, index: bool = False, key: str | None = None):
    csv = df.to_csv(index=index).encode("utf-8")
    st.download_button(label, data=csv, file_name=filename, mime="text/csv", key=key)


# --------------------------------------------------------------------------
# Sidebar: data input
# --------------------------------------------------------------------------

st.sidebar.title("\U0001F9EC SatSub")
st.sidebar.caption("MSA statistics & nucleotide substitution saturation")

st.sidebar.header("1. Nucleotide alignment")
nt_source = st.sidebar.radio(
    "Source",
    ["Upload a file", "Use bundled example (Vertebrate COI, 8 taxa)"],
    key="nt_source",
    label_visibility="collapsed",
)

nt_loaded = None
if nt_source == "Use bundled example (Vertebrate COI, 8 taxa)":
    nt_loaded = read_alignment(example_alignment_path())
else:
    nt_file = st.sidebar.file_uploader(
        "FASTA / Clustal / Phylip / Nexus / Stockholm",
        type=["fasta", "fa", "fas", "fna", "aln", "clustal", "phy", "phylip", "nex", "nexus", "sto", "stockholm", "txt"],
        key="nt_upload",
    )
    if nt_file is not None:
        try:
            nt_loaded = read_alignment(nt_file)
            if nt_loaded.kind != "nucleotide":
                st.sidebar.warning(
                    "This file looks like amino-acid data, not nucleotides. "
                    "Upload it below as the amino-acid alignment instead."
                )
        except AlignmentReadError as exc:
            st.sidebar.error(str(exc))

frame_start = 1
if nt_loaded is not None and nt_loaded.kind == "nucleotide":
    st.sidebar.header("2. Reading frame")
    frame_start = st.sidebar.selectbox(
        "Codon position 1 starts at alignment column",
        options=[1, 2, 3],
        index=0,
        help=(
            "Used for codon-position GC content, for splitting s/v-vs-distance "
            "plots by codon position, and for translation. Codon position is "
            "assigned by column index only; indels that break the reading "
            "frame are not re-detected."
        ),
    )

st.sidebar.header("3. Amino-acid alignment")
aa_default_index = 0 if (nt_loaded is not None and nt_loaded.kind == "nucleotide") else 1
aa_source = st.sidebar.radio(
    "Source",
    ["Translate the nucleotide alignment", "Upload a separate amino-acid alignment", "Skip"],
    index=aa_default_index,
    key="aa_source",
    label_visibility="collapsed",
)

aa_pairs = None
aa_warnings: list[str] = []
if aa_source == "Translate the nucleotide alignment":
    if nt_loaded is None or nt_loaded.kind != "nucleotide":
        st.sidebar.info("Load a nucleotide alignment first.")
    else:
        table_ids = list(GENETIC_CODE_CHOICES)
        table_id = st.sidebar.selectbox(
            "Genetic code table",
            options=table_ids,
            index=table_ids.index(1),
            format_func=lambda i: f"{i} – {GENETIC_CODE_CHOICES[i]}",
            help="Vertebrate/invertebrate mitochondrial genes (e.g. COI barcoding) typically need table 2 or 5.",
        )
        try:
            aa_pairs, aa_warnings = _cached_translate(
                alignment_to_pairs(nt_loaded.alignment), frame_start, table_id
            )
        except ValueError as exc:
            st.sidebar.error(str(exc))
elif aa_source == "Upload a separate amino-acid alignment":
    aa_file = st.sidebar.file_uploader(
        "FASTA / Clustal / Phylip / Nexus / Stockholm",
        type=["fasta", "fa", "fas", "faa", "aln", "clustal", "phy", "phylip", "nex", "nexus", "sto", "stockholm", "txt"],
        key="aa_upload",
    )
    if aa_file is not None:
        try:
            aa_loaded = read_alignment(aa_file)
            if aa_loaded.kind != "amino_acid":
                st.sidebar.warning("This file looks like nucleotide data, not amino acids.")
            aa_pairs = alignment_to_pairs(aa_loaded.alignment)
            aa_warnings = aa_loaded.warnings
        except AlignmentReadError as exc:
            st.sidebar.error(str(exc))

st.sidebar.divider()
st.sidebar.caption("SatSub v0.1 · Biopython + NumPy/SciPy + Streamlit/Plotly")


# --------------------------------------------------------------------------
# Main area
# --------------------------------------------------------------------------

st.title("SatSub: MSA statistics & substitution saturation")

tab_overview, tab_nt, tab_aa, tab_sat, tab_about = st.tabs(
    ["Overview", "Nucleotide statistics", "Amino-acid statistics", "Substitution saturation", "About / methods"]
)

with tab_overview:
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Nucleotide alignment")
        if nt_loaded is None:
            st.info("Load a nucleotide alignment from the sidebar to get started.")
        else:
            st.write(f"**Source:** {nt_loaded.source_name}")
            st.write(f"**Format detected:** {nt_loaded.format}")
            st.write(f"**Sequences:** {nt_loaded.n_sequences}  ·  **Length:** {nt_loaded.length} bp")
            if nt_loaded.had_rna:
                st.caption("Input used the RNA alphabet (U); normalized to T internally.")
            for w in nt_loaded.warnings:
                st.warning(w)
    with col2:
        st.subheader("Amino-acid alignment")
        if aa_pairs is None:
            st.info("No amino-acid alignment loaded (translate or upload one from the sidebar).")
        else:
            aa_aln_preview = alignment_from_pairs(aa_pairs)
            st.write(f"**Sequences:** {len(aa_aln_preview)}  ·  **Length:** {aa_aln_preview.get_alignment_length()} aa")
            for w in aa_warnings:
                st.warning(w)

with tab_nt:
    if nt_loaded is None or nt_loaded.kind != "nucleotide":
        st.info("Load a nucleotide alignment from the sidebar to see statistics here.")
    else:
        summary = _cached_nt_summary(alignment_to_pairs(nt_loaded.alignment), frame_start)
        o = summary.overall

        st.subheader("Summary")
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Sequences", o["n_sequences"])
        c2.metric("Alignment length", o["alignment_length"])
        c3.metric("Overall GC%", f"{o['overall_GC_percent']:.1f}")
        c4.metric("Overall gaps%", f"{o['overall_gap_percent']:.1f}")
        c5, c6, c7, c8 = st.columns(4)
        c5.metric("Mean pairwise identity%", f"{o['mean_pairwise_identity_percent']:.1f}")
        c6.metric("Mean pairwise p-distance", f"{o['mean_pairwise_p_distance']:.4f}")
        c7.metric("Aggregate ts/tv ratio", f"{o['aggregate_ts_tv_ratio']:.2f}")
        c8.metric("Ambiguous/other%", f"{o['ambiguous_or_other_percent']:.2f}")

        st.subheader("Site classification")
        c1, c2, c3, c4, c5 = st.columns(5)
        c1.metric("Conserved", o["conserved_sites"])
        c2.metric("Variable", o["variable_sites"])
        c3.metric("Parsimony-informative", o["parsimony_informative_sites"])
        c4.metric("Singleton", o["singleton_sites"])
        c5.metric("Fully gapped", o["fully_gapped_sites"])

        if "GC1_percent" in o:
            st.subheader("Codon-position GC content")
            c1, c2, c3, c4 = st.columns(4)
            c1.metric("GC1%", f"{o['GC1_percent']:.1f}")
            c2.metric("GC2%", f"{o['GC2_percent']:.1f}")
            c3.metric("GC3%", f"{o['GC3_percent']:.1f}")
            c4.metric("GC12 mean%", f"{o['GC12_mean_percent']:.1f}")

        col1, col2 = st.columns(2)
        with col1:
            st.plotly_chart(
                categorical_bar(
                    ["A", "C", "G", "T"],
                    [summary.per_sequence[b].sum() for b in "ACGT"],
                    "Overall base composition",
                    "Count",
                ),
                width="stretch",
            )
        with col2:
            st.plotly_chart(identity_heatmap(summary.identity_matrix), width="stretch")

        st.subheader("Per-sequence statistics")
        st.dataframe(summary.per_sequence, width="stretch")
        _df_download_button(summary.per_sequence, "Download per-sequence CSV", "satsub_nt_per_sequence.csv", key="nt_perseq_dl")
        _df_download_button(summary.overall_table(), "Download overall-summary CSV", "satsub_nt_overall.csv", key="nt_overall_dl")
        _df_download_button(
            summary.identity_matrix, "Download identity-matrix CSV", "satsub_nt_identity.csv", index=True, key="nt_ident_dl"
        )

with tab_aa:
    if aa_pairs is None:
        st.info("Translate the nucleotide alignment or upload an amino-acid alignment from the sidebar.")
    else:
        aa_summary = _cached_aa_summary(aa_pairs)
        o = aa_summary.overall

        st.subheader("Summary")
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Sequences", o["n_sequences"])
        c2.metric("Alignment length", o["alignment_length"])
        c3.metric("Overall gaps%", f"{o['overall_gap_percent']:.1f}")
        c4.metric("Ambiguous/other%", f"{o['ambiguous_or_other_percent']:.2f}")
        c5, c6 = st.columns(2)
        c5.metric("Mean pairwise identity%", f"{o['mean_pairwise_identity_percent']:.1f}")
        c6.metric("Mean pairwise similarity% (BLOSUM62>0)", f"{o['mean_pairwise_similarity_percent']:.1f}")

        st.subheader("Site classification")
        c1, c2, c3, c4, c5 = st.columns(5)
        c1.metric("Conserved", o["conserved_sites"])
        c2.metric("Variable", o["variable_sites"])
        c3.metric("Parsimony-informative", o["parsimony_informative_sites"])
        c4.metric("Singleton", o["singleton_sites"])
        c5.metric("Fully gapped", o["fully_gapped_sites"])

        st.subheader("Amino-acid composition")
        aa_letters = list("ACDEFGHIKLMNPQRSTVWY")
        st.plotly_chart(
            monochrome_bar(
                aa_letters,
                [
                    sum(str(rec.seq).count(a) for rec in alignment_from_pairs(aa_pairs))
                    for a in aa_letters
                ],
                "Overall amino-acid composition",
                "Count",
            ),
            width="stretch",
        )

        st.subheader("Physicochemical composition (mean of per-sequence %)")
        physico_cats = ["hydrophobic_percent", "polar_uncharged_percent", "acidic_percent", "basic_percent"]
        st.plotly_chart(
            categorical_bar(
                ["Hydrophobic", "Polar uncharged", "Acidic", "Basic"],
                [aa_summary.per_sequence[c].mean() for c in physico_cats],
                "Physicochemical group composition",
                "Mean % of classified residues",
            ),
            width="stretch",
        )

        col1, col2 = st.columns(2)
        with col1:
            st.plotly_chart(identity_heatmap(aa_summary.identity_matrix, "Pairwise identity (%)"), width="stretch")
        with col2:
            st.plotly_chart(
                identity_heatmap(aa_summary.similarity_matrix, "Pairwise similarity (%, BLOSUM62 > 0)"),
                width="stretch",
            )

        st.subheader("Per-sequence statistics")
        st.dataframe(aa_summary.per_sequence, width="stretch")
        _df_download_button(aa_summary.per_sequence, "Download per-sequence CSV", "satsub_aa_per_sequence.csv", key="aa_perseq_dl")
        _df_download_button(aa_summary.overall_table(), "Download overall-summary CSV", "satsub_aa_overall.csv", key="aa_overall_dl")

with tab_sat:
    if nt_loaded is None or nt_loaded.kind != "nucleotide":
        st.info("Load a nucleotide alignment from the sidebar to run the saturation analysis.")
    elif nt_loaded.n_sequences < 2:
        st.warning("Need at least 2 sequences to compute pairwise distances.")
    else:
        st.markdown(
            "For every sequence pair, observed transition (**s**) and transversion (**v**) "
            "proportions are plotted against a model-corrected genetic distance. "
            "If points flatten out (plateau) as distance increases instead of tracking it "
            "linearly, that pattern signals **substitution saturation** at those sites: "
            "multiple hits have started to erase the phylogenetic signal, and that distance "
            "correction/site class should be interpreted cautiously in tree-building."
        )
        pos_tabs = st.tabs([POSITION_LABELS[k] for k in ("all", "1", "2", "3")])
        for label_key, pos_tab in zip(("all", "1", "2", "3"), pos_tabs):
            with pos_tab:
                try:
                    result = _cached_saturation(alignment_to_pairs(nt_loaded.alignment), label_key, frame_start)
                except ValueError as exc:
                    st.error(str(exc))
                    continue
                df = result.table
                c1, c2, c3, c4 = st.columns(4)
                c1.metric("Sites used", result.n_sites_used)
                c2.metric("Pairs", len(df))
                c3.metric("Mean p-distance", f"{df['p_distance'].mean():.4f}")
                c4.metric("Mean ts/tv ratio", f"{df['ts_tv_ratio'].replace([float('inf')], float('nan')).mean():.2f}")

                grid = [("jc69", "k80"), ("tn93", "gtr")]
                for row_models in grid:
                    gcol1, gcol2 = st.columns(2)
                    for model_key, gcol in zip(row_models, (gcol1, gcol2)):
                        with gcol:
                            fig, n_dropped = saturation_scatter(
                                df, model_key, MODEL_LABELS[model_key], MODEL_LABELS[model_key]
                            )
                            st.plotly_chart(fig, width="stretch")
                            if n_dropped:
                                st.caption(
                                    f"{n_dropped} of {len(df)} pairs omitted: distance undefined "
                                    "under this model (fully saturated at these sites)."
                                )

                with st.expander("Raw pairwise data"):
                    st.dataframe(df, width="stretch")
                    _df_download_button(
                        df,
                        f"Download CSV ({POSITION_LABELS[label_key]})",
                        f"satsub_saturation_{label_key}.csv",
                        key=f"sat_dl_{label_key}",
                    )

with tab_about:
    st.markdown(
        """
### What SatSub computes

**Alignment statistics** (Nucleotide / Amino-acid tabs): sequence counts and
lengths, gap and ambiguous-character percentages, base/amino-acid
composition, GC content (including per-codon-position GC1/GC2/GC3 for
coding nucleotide alignments), conserved / variable / parsimony-informative
/ singleton site counts, and pairwise identity (amino acids also get a
BLOSUM62-based similarity score, and a simple hydrophobic / polar / acidic
/ basic composition breakdown).

**Substitution saturation** (Saturation tab): for every sequence pair, sites
are classified as identical, a transition (purine A↔G or pyrimidine C↔T),
or a transversion, optionally restricted to a single codon position. Four
distance corrections are then computed from those counts:

- **JC69** (Jukes & Cantor, 1969) – equal base frequencies, equal
  substitution rates.
- **K80** (Kimura, 1980) – separate transition/transversion rates, equal
  base frequencies.
- **TN93** (Tamura & Nei, 1993) – unequal base frequencies, separate
  purine-transition, pyrimidine-transition and transversion rates.
- **GTR** (general time-reversible) – unequal base frequencies and all six
  substitution types rated independently.

Plotting each pair's observed **s** and **v** proportions against the
model-corrected distance reproduces the classic saturation plot (e.g. Xia
2003, Philippe 1994): when substitutions accumulate faster than the
corrected distance can track (multiple/back substitutions at the same
site), the points plateau instead of increasing linearly, which flags that
site class / model as unreliable for deep divergences.

#### A note on the GTR distance

Unlike JC69, K80 and TN93, GTR has no closed-form pairwise distance
formula. SatSub estimates the six exchangeability rates once per alignment
(or per codon-position subset) from substitution counts pooled across every
sequence pair, then computes each pair's distance as the maximum-likelihood
branch length under that fixed rate matrix. This "estimate once, then
profile per pair" strategy is the standard way to get pairwise GTR-like
distances without fitting a full tree, but it is an aggregate/composite
approximation rather than a joint phylogenetic ML fit — treat it as
exploratory, the same way the JC69/K80/TN93 saturation plots are meant to
be exploratory rather than a final phylogenetic analysis.

#### Codon positions

Codon position is assigned purely from alignment-column index and the
reading-frame start chosen in the sidebar. It does not re-derive frames
around indels, so alignments with frameshifting gaps (length not a
multiple of 3) will have codon positions shift for all downstream columns;
use a codon-aware alignment for such data.

#### References

Jukes TH & Cantor CR (1969). *Evolution of Protein Molecules.* In Mammalian
Protein Metabolism, 21–32.
Kimura M (1980). A simple method for estimating evolutionary rates of base
substitutions. *J Mol Evol* 16:111–120.
Tamura K & Nei M (1993). Estimation of the number of nucleotide
substitutions in the control region of mitochondrial DNA in humans and
chimpanzees. *Mol Biol Evol* 10:512–526.
Tavaré S (1986). Some probabilistic and statistical problems in the
analysis of DNA sequences. *Lect Math Life Sci* 17:57–86.
        """
    )
