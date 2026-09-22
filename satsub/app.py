"""SatSub Streamlit GUI.

Run with:  streamlit run satsub/app.py
or, once installed:  satsub-gui
"""

from __future__ import annotations

import pandas as pd
import streamlit as st

from satsub.branding import LAB_NAME, LAB_URL, LOGO_PATH
from satsub.composition import (
    codon_usage_table,
    directional_pair_frequencies,
    nucleotide_composition_table,
)
from satsub.io import (
    GENETIC_CODE_CHOICES,
    AlignmentReadError,
    alignment_from_pairs,
    alignment_to_pairs,
    example_alignment_path,
    read_alignment,
    translate_alignment,
)
from satsub.plotting import (
    categorical_bar,
    directional_pair_heatmap,
    grouped_bar_by_base,
    identity_heatmap,
    monochrome_bar,
    rscu_diverging_bar,
    saturation_scatter,
)
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


@st.cache_data(show_spinner="Computing nucleotide composition...")
def _cached_nt_composition(pairs, frame_start):
    aln = alignment_from_pairs(pairs)
    return nucleotide_composition_table(aln, frame_start=frame_start)


@st.cache_data(show_spinner="Computing codon usage...")
def _cached_codon_usage(pairs, frame_start, table_id):
    aln = alignment_from_pairs(pairs)
    return codon_usage_table(aln, frame_start=frame_start, table_id=table_id)


@st.cache_data(show_spinner="Computing directional pair frequencies...")
def _cached_directional_pairs(pairs, frame_start):
    aln = alignment_from_pairs(pairs)
    return directional_pair_frequencies(aln, frame_start=frame_start)


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

def _image_data_uri(path) -> str | None:
    try:
        data = path.read_bytes()
    except OSError:
        return None
    import base64

    return "data:image/png;base64," + base64.b64encode(data).decode()


_logo_uri = _image_data_uri(LOGO_PATH)
if _logo_uri is not None:
    st.sidebar.markdown(
        f'<a href="{LAB_URL}" target="_blank" rel="noopener">'
        f'<img src="{_logo_uri}" alt="{LAB_NAME}" style="width:100%; max-width:260px; margin-bottom:0.5rem;">'
        f"</a>",
        unsafe_allow_html=True,
    )

st.sidebar.title("\U0001F9EC SatSub")
st.sidebar.caption("MSA statistics & nucleotide substitution saturation")
st.sidebar.markdown(f"A tool from the [{LAB_NAME}]({LAB_URL})")

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
table_id = 1
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

    st.sidebar.header("3. Genetic code table")
    table_ids = list(GENETIC_CODE_CHOICES)
    table_id = st.sidebar.selectbox(
        "Used for translation, codon usage and RSCU",
        options=table_ids,
        index=table_ids.index(1),
        format_func=lambda i: f"{i} – {GENETIC_CODE_CHOICES[i]}",
        help="Vertebrate/invertebrate mitochondrial genes (e.g. COI barcoding) typically need table 2 or 5.",
        key="genetic_code_table_select",
    )
    st.sidebar.caption(
        "This single choice is shared by amino-acid translation (below) and by "
        "the Composition & codon usage tab, so the two never disagree."
    )

st.sidebar.header("4. Amino-acid alignment")
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
st.sidebar.caption(f"[{LAB_URL.removeprefix('https://')}]({LAB_URL})")


# --------------------------------------------------------------------------
# Main area
# --------------------------------------------------------------------------

st.title("SatSub: MSA statistics & substitution saturation")

tab_overview, tab_nt, tab_comp, tab_aa, tab_sat, tab_about = st.tabs(
    [
        "Overview",
        "Nucleotide statistics",
        "Composition & codon usage",
        "Amino-acid statistics",
        "Substitution saturation",
        "About / methods",
    ]
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

with tab_comp:
    if nt_loaded is None or nt_loaded.kind != "nucleotide":
        st.info("Load a nucleotide alignment from the sidebar to see composition and codon usage here.")
    else:
        nt_pairs = alignment_to_pairs(nt_loaded.alignment)

        st.subheader("Nucleotide frequencies")
        st.caption(
            "Per-sequence T/C/A/G composition, overall and at each codon position. "
            "The chart shows the mean across sequences (the 'Avg.' table row)."
        )
        nc = _cached_nt_composition(nt_pairs, frame_start)
        avg_row = nc[nc["sequence"] == "Avg."].iloc[0]
        if "T1_percent" in nc.columns:
            x_cats = ["All", "1st", "2nd", "3rd"]
            series = {
                b: [
                    avg_row[f"{b}_percent"],
                    avg_row[f"{b}1_percent"],
                    avg_row[f"{b}2_percent"],
                    avg_row[f"{b}3_percent"],
                ]
                for b in "ACGT"
            }
        else:
            x_cats = ["All"]
            series = {b: [avg_row[f"{b}_percent"]] for b in "ACGT"}
        st.plotly_chart(
            grouped_bar_by_base(x_cats, series, "Mean nucleotide composition by codon position", "%"),
            width="stretch",
        )
        with st.expander("Per-sequence nucleotide frequency table"):
            st.dataframe(nc, width="stretch")
            _df_download_button(nc, "Download CSV", "satsub_nucleotide_composition.csv", key="nt_comp_dl")

        st.divider()
        st.subheader("Codon usage")
        st.caption(
            f"Using genetic code table {table_id} – {GENETIC_CODE_CHOICES[table_id]}, "
            "as selected in the sidebar (shared with amino-acid translation). "
            "Change it there to recompute this section under a different code."
        )
        try:
            cu = _cached_codon_usage(nt_pairs, frame_start, table_id)
            c1, c2 = st.columns(2)
            c1.metric("Mean codons per sequence", f"{cu.avg_codons_per_sequence:.1f}")
            c2.metric("Genetic code", cu.table_name)
            st.plotly_chart(
                monochrome_bar(
                    [f"{r.codon} ({r.amino_acid})" for r in cu.table.itertuples()],
                    list(cu.table["count_avg"]),
                    "Mean codon usage (average count per sequence)",
                    "Mean count",
                ).update_xaxes(tickangle=-90, tickfont=dict(size=9)),
                width="stretch",
            )
            st.plotly_chart(rscu_diverging_bar(cu.table), width="stretch")
            st.caption(
                "RSCU (relative synonymous codon usage) compares each codon's usage to the average "
                "usage of its synonymous family under the selected genetic code; 1.0 (dotted line) "
                "is the value expected if all synonyms were used equally."
            )
            with st.expander("Codon usage table (all 64 codons)"):
                st.dataframe(cu.table, width="stretch")
                _df_download_button(cu.table, "Download CSV", "satsub_codon_usage.csv", key="codon_usage_dl")
        except ValueError as exc:
            st.error(str(exc))

        st.divider()
        st.subheader("Directional base-pair frequencies")
        st.caption(
            "For every sequence pair and site, the aligned base pair is classified as identical (ii), "
            "a transition (si) or a transversion (sv), then averaged over every pairwise comparison. "
            "Row = base in the earlier-listed sequence of each pair; column = base in the later one "
            "(an arbitrary but fixed direction, since pooling over unordered pairs has no natural order)."
        )
        dp = _cached_directional_pairs(nt_pairs, frame_start)
        summary_display = dp.summary.copy()
        summary_display["subset"] = summary_display["subset"].map(POSITION_LABELS)
        st.dataframe(summary_display, width="stretch", hide_index=True)
        _df_download_button(dp.summary, "Download summary CSV", "satsub_directional_pairs_summary.csv", key="dp_summary_dl")

        pair_tabs = st.tabs([POSITION_LABELS[k] for k in ("all", "1", "2", "3") if k in dp.matrices])
        for key, ptab in zip([k for k in ("all", "1", "2", "3") if k in dp.matrices], pair_tabs):
            with ptab:
                st.plotly_chart(
                    directional_pair_heatmap(dp.matrices[key], f"Directional pairs – {POSITION_LABELS[key]}"),
                    width="stretch",
                )
                _df_download_button(
                    dp.matrices[key],
                    "Download matrix CSV",
                    f"satsub_directional_pairs_{key}.csv",
                    index=True,
                    key=f"dp_matrix_dl_{key}",
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
    if _logo_uri is not None:
        st.markdown(
            f'<a href="{LAB_URL}" target="_blank" rel="noopener">'
            f'<img src="{_logo_uri}" alt="{LAB_NAME}" style="max-width:220px;">'
            f"</a>",
            unsafe_allow_html=True,
        )
    st.caption(f"SatSub is developed by the [{LAB_NAME}]({LAB_URL}) ({LAB_URL.removeprefix('https://')}).")
    st.caption(
        "Licensed under the [MIT License](https://opensource.org/license/mit) "
        "– free to use, modify and redistribute, including commercially, "
        "with the copyright notice preserved."
    )

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

**Composition & codon usage** (Composition & codon usage tab, nucleotide
alignments only):

- *Nucleotide frequencies* – T/C/A/G percentages per sequence, overall and
  at each codon position, with an "Avg." row/chart giving the simple mean
  across sequences.
- *Codon usage* – for each of the 64 codons, the mean raw count per
  sequence (averaged across taxa, not a per-mille rate) and its **relative
  synonymous codon usage (RSCU)**: the codon's mean count divided by the
  average count of its synonymous family under the selected genetic code
  (Sharp, Tuohy & Mosurski, 1986). RSCU = 1 means the codon is used exactly
  as often as expected if every synonym were used equally; the RSCU chart's
  diverging blue/red color is centered on that value.
- *Directional base-pair frequencies* – for every sequence pair and site,
  the aligned base pair is classified as identical (**ii**), a transition
  (**si**) or a transversion (**sv**), then averaged across every pairwise
  comparison (giving **R = si/sv**), plus the full averaged 4×4 substitution
  matrix. Because pooling over many unordered sequence pairs has no natural
  "first vs. second" direction, SatSub fixes one convention throughout: the
  row is the base in whichever sequence is listed earlier in the alignment,
  the column is the base in the later one. This is an arbitrary but
  internally consistent labeling, not a claim about ancestry or direction
  of change.

All three are computed for the whole alignment and, where a reading frame
is available, separately for each codon position.

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
site class / model as unreliable for deep divergences. Each series also
gets a **LOWESS** (locally weighted) smoothed trend line rather than a
single straight fit, so a plateau shows up as a visible bend instead of
being averaged away.

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

#### Genetic code table

The sidebar has a single "Genetic code table" choice, used both to
translate the nucleotide alignment to amino acids and to compute codon
usage/RSCU in the Composition & codon usage tab. Changing it updates both
consistently, so translation and codon usage can never disagree about which
code was used.

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
Sharp PM, Tuohy TM & Mosurski KR (1986). Codon usage in yeast: cluster
analysis clearly differentiates highly and lowly expressed genes. *Nucleic
Acids Res* 14:5125–5143.
Cleveland WS (1979). Robust locally weighted regression and smoothing
scatterplots. *J Am Stat Assoc* 74:829–836.
        """
    )
