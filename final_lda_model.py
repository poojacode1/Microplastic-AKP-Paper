import marimo

__generated_with = "0.23.10"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    mo.md(
        """
        # Microbiome LDA Topic Modeling

        """
    )
    return


@app.cell
def _():
    import warnings

    warnings.filterwarnings("ignore")

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib.colors as mcolors
    import seaborn as sns

    from gensim import corpora, models
    from gensim.models import CoherenceModel
    from scipy.cluster.hierarchy import dendrogram, linkage

    return (
        CoherenceModel,
        corpora,
        dendrogram,
        linkage,
        mcolors,
        models,
        mpatches,
        np,
        pd,
        plt,
        sns,
    )


@app.cell
def _():

    META_PATH = "meta-edited.csv"
    OTU_PATH = "otu_table.csv"
    TAX_PATH = "taxonomy_assignments.csv"

    MIN_TOTAL_READS = 10   # drop ASVs with < 10 reads across the whole dataset
    MIN_READS = 100       

    TOPIC_MIN = 2
    TOPIC_MAX = 13
    RANDOM_STATE = 42
    LDA_PASSES = 10
    return (
        LDA_PASSES,
        META_PATH,
        MIN_READS,
        MIN_TOTAL_READS,
        OTU_PATH,
        RANDOM_STATE,
        TAX_PATH,
        TOPIC_MAX,
        TOPIC_MIN,
    )


@app.cell
def _(META_PATH, pd):
    #load metadata
    metadata = pd.read_csv(META_PATH, sep=";", encoding="utf-8-sig")

    print(f"Metadata: {metadata.shape[0]} samples")
    print(f"Columns: {', '.join(metadata.columns)}")
    metadata.head()
    return (metadata,)


@app.cell
def _(OTU_PATH, pd):
    #load asv table
    otu_raw = pd.read_csv(OTU_PATH, index_col=0)

    print(f"OTU table: {otu_raw.shape[0]} samples x {otu_raw.shape[1]} OTUs")
    otu_raw.iloc[:5, :5]
    return (otu_raw,)


@app.cell
def _(TAX_PATH, pd):
    #load taxonomy
    taxonomy = pd.read_csv(TAX_PATH, index_col=0)

    print(f"Taxonomy: {taxonomy.shape[0]} taxa")
    print(f"Columns: {', '.join(taxonomy.columns)}")
    taxonomy.head()
    return (taxonomy,)


@app.cell
def _(metadata):
    #Build the list of samples to discard 
    # Daphnia env. samples, kit controls, sterilised water, uncertain plastic
    # assignment, and all bacterioplankton.
    samples_to_remove = []

    if "Category" in metadata.columns:
        _daphnia_env = metadata.loc[
            metadata["Category"] == "Daphnia:env_sample", "sample_name"
        ].tolist()
        if _daphnia_env:
            samples_to_remove.extend(_daphnia_env)
            print(f"- Daphnia:env_sample: {len(_daphnia_env)}")

        _kit_controls = metadata.loc[
            metadata["Category"] == "Kit-control", "sample_name"
        ].tolist()
        if _kit_controls:
            samples_to_remove.extend(_kit_controls)
            print(f"- Kit-control: {len(_kit_controls)}")

        if "Sample _info" in metadata.columns:
            _steril = metadata.loc[
                (metadata["Category"] == "Bacterioplankton: Post-exposure")
                & metadata["Sample _info"].str.contains("steril", case=False, na=False),
                "sample_name",
            ].tolist()
            if _steril:
                samples_to_remove.extend(_steril)
                print(f"- Sterilised water: {len(_steril)}")

        if "Plastic" in metadata.columns:
            _not_sure = metadata.loc[
                metadata["Category"].str.contains("Daphnia", case=False, na=False)
                & metadata["Plastic"].str.contains("not sure", case=False, na=False),
                "sample_name",
            ].tolist()
            if _not_sure:
                samples_to_remove.extend(_not_sure)
                print(f"- Daphnia 'not sure' plastic: {len(_not_sure)}")

        _bacterio = metadata.loc[
            metadata["Category"].str.contains("Bacterioplankton", case=False, na=False),
            "sample_name",
        ].tolist()
        if _bacterio:
            samples_to_remove.extend(_bacterio)
            print(f"- Bacterioplankton: {len(_bacterio)}")
    else:
        print("No 'Category' column found - nothing removed.")

    samples_to_remove = sorted(set(samples_to_remove))
    print(f"\nTotal unique samples flagged for removal: {len(samples_to_remove)}")
    return (samples_to_remove,)


@app.cell
def _(metadata, otu_raw, samples_to_remove):
    # Apply the control/ sample filter 
    otu_no_ctrl = otu_raw[~otu_raw.index.isin(samples_to_remove)].copy()
    meta_no_ctrl = metadata[~metadata["sample_name"].isin(samples_to_remove)].copy()

    print(f"Remaining samples: {otu_no_ctrl.shape[0]}")
    return meta_no_ctrl, otu_no_ctrl


@app.cell
def _(MIN_TOTAL_READS, otu_no_ctrl):
    #drop ASVs with < 10 reads total
    _nz = otu_no_ctrl.loc[:, otu_no_ctrl.sum(axis=0) > 0]
    print(f"After removing zero-count ASVs: {_nz.shape[1]} ASVs")

    _keep = _nz.sum(axis=0) >= MIN_TOTAL_READS
    otu_prevalent = _nz.loc[:, _keep]

    print(f"After frequency filter (>= {MIN_TOTAL_READS} reads total): "
          f"{otu_prevalent.shape[1]} ASVs")
    print(f"Reads retained: "
          f"{100 * otu_prevalent.values.sum() / _nz.values.sum():.2f}%")
    return (otu_prevalent,)


@app.cell
def _(MIN_READS, meta_no_ctrl, otu_prevalent):
    #drop low depth samples
    _totals = otu_prevalent.sum(axis=1)
    _low = _totals[_totals < MIN_READS].index.tolist()

    if _low:
        print(f"Removing {len(_low)} samples with < {MIN_READS} reads")
        otu_deep = otu_prevalent[~otu_prevalent.index.isin(_low)]
        meta_deep = meta_no_ctrl[~meta_no_ctrl["sample_name"].isin(_low)]
    else:
        print("No low-count samples to remove")
        otu_deep = otu_prevalent
        meta_deep = meta_no_ctrl

    print(f"Dataset: {otu_deep.shape[0]} samples x {otu_deep.shape[1]} OTUs")
    return meta_deep, otu_deep


@app.cell
def _(meta_deep, otu_deep):
    #Align OTU table and metadata on a common sample set 
    common_samples = sorted(set(otu_deep.index) & set(meta_deep["sample_name"]))

    otu_final = otu_deep.loc[common_samples]
    meta_final = (
        meta_deep.set_index("sample_name").loc[common_samples].reset_index()
    )

    print(f"Aligned on {len(common_samples)} samples")
    return meta_final, otu_final


@app.cell
def _(otu_final, pd, taxonomy):
    #Readable taxonomy names 
    def create_otu_name(row):
        """Build a readable label from the finest available rank."""
        if pd.notna(row.get("Genus")):
            _genus = row["Genus"]
            if pd.notna(row.get("Species")):
                return f"{_genus}_{row['Species']}"
            return f"{_genus}_sp"
        for _rank in ("Family", "Order", "Class", "Phylum"):
            if pd.notna(row.get(_rank)):
                return f"{row[_rank]}_uncl"
        return "Unclassified"

    taxonomy_filtered = taxonomy.loc[taxonomy.index.isin(otu_final.columns)].copy()
    taxonomy_filtered["readable_name"] = taxonomy_filtered.apply(create_otu_name, axis=1)

    # Disambiguate duplicated labels
    _counts = taxonomy_filtered["readable_name"].value_counts()
    for _name in _counts[_counts > 1].index:
        _idxs = taxonomy_filtered.index[taxonomy_filtered["readable_name"] == _name]
        for _i, _idx in enumerate(_idxs, start=1):
            taxonomy_filtered.loc[_idx, "readable_name"] = f"{_name}_{_i}"

    otu_name_mapping = taxonomy_filtered["readable_name"].to_dict()
    print(f"Created readable names for {len(otu_name_mapping)} OTUs")
    return otu_name_mapping, taxonomy_filtered


@app.cell
def _(corpora, otu_final, otu_name_mapping):
    # Documents, dictionary, corpus 
    documents = []
    for _sample_id, _row in otu_final.iterrows():
        _doc = []
        for _otu_seq, _count in _row.items():
            if _count > 0:
                _name = otu_name_mapping.get(_otu_seq, str(_otu_seq)[:20])
                _doc.extend([_name] * int(round(_count)))
        documents.append(_doc)

    dictionary = corpora.Dictionary(documents)
    corpus = [dictionary.doc2bow(_d) for _d in documents]

    print(f"Corpus: {len(corpus)} documents")
    print(f"Dictionary: {len(dictionary)} unique taxa")
    return corpus, dictionary, documents


@app.cell
def _(
    CoherenceModel,
    LDA_PASSES,
    RANDOM_STATE,
    TOPIC_MAX,
    TOPIC_MIN,
    corpus,
    dictionary,
    documents,
    models,
):
    #finding optimal topics
    def compute_coherence_values(dictionary, corpus, texts, start, limit, step=1):
        _models, _coherences, _perplexities = [], [], []

        for _k in range(start, limit, step):
            print(f"   Fitting {_k} topics...", end="\r")
            _model = models.LdaModel(
                corpus=corpus,
                id2word=dictionary,
                num_topics=_k,
                random_state=RANDOM_STATE,
                passes=LDA_PASSES,
                alpha="auto",
                per_word_topics=True,
            )
            _models.append(_model)

            _cm = CoherenceModel(
                model=_model, texts=texts, dictionary=dictionary, coherence="c_v"
            )
            _coherences.append(_cm.get_coherence())
            _perplexities.append(_model.log_perplexity(corpus))

        print()
        return _models, _coherences, _perplexities

    model_list, coherence_values, perplexity_values = compute_coherence_values(
        dictionary=dictionary,
        corpus=corpus,
        texts=documents,
        start=TOPIC_MIN,
        limit=TOPIC_MAX,
        step=1,
    )
    return coherence_values, model_list, perplexity_values


@app.cell
def _(TOPIC_MIN, coherence_values, model_list, np):
    # Pick the best model 
    optimal_num_topics = int(np.argmax(coherence_values)) + TOPIC_MIN
    lda_model = model_list[optimal_num_topics - TOPIC_MIN]
    best_coherence = coherence_values[optimal_num_topics - TOPIC_MIN]

    print(f"Optimal number of topics: {optimal_num_topics}")
    print(f"Best coherence score: {best_coherence:.4f}")
    return best_coherence, lda_model, optimal_num_topics


@app.cell
def _(
    TOPIC_MAX,
    TOPIC_MIN,
    coherence_values,
    optimal_num_topics,
    perplexity_values,
    plt,
):
    #  Coherence / perplexity curves 
    _fig, (_ax1, _ax2) = plt.subplots(1, 2, figsize=(14, 5))
    _x = range(TOPIC_MIN, TOPIC_MAX)

    _ax1.plot(_x, coherence_values, marker="o", linewidth=2, markersize=8)
    _ax1.axvline(
        optimal_num_topics, color="r", linestyle="--",
        label=f"Optimal: {optimal_num_topics}",
    )
    _ax1.set_xlabel("Number of Topics", fontsize=12)
    _ax1.set_ylabel("Coherence Score", fontsize=12)
    _ax1.set_title("Model Coherence", fontsize=14, fontweight="bold")
    _ax1.legend()
    _ax1.grid(True, alpha=0.3)

    _ax2.plot(
        _x, perplexity_values, marker="s", linewidth=2, markersize=8, color="orange"
    )
    _ax2.axvline(
        optimal_num_topics, color="r", linestyle="--",
        label=f"Optimal: {optimal_num_topics}",
    )
    _ax2.set_xlabel("Number of Topics", fontsize=12)
    _ax2.set_ylabel("Log Perplexity", fontsize=12)
    _ax2.set_title("Model Perplexity", fontsize=14, fontweight="bold")
    _ax2.legend()
    _ax2.grid(True, alpha=0.3)

    _fig.tight_layout()
    _fig.savefig("model_selection_metrics.png", dpi=300, bbox_inches="tight")
    _fig
    return


@app.cell
def _(corpus, lda_model, optimal_num_topics, otu_final, pd):
    #Per-sample topic distributions 
    _distributions, _dominant = [], []

    for _bow in corpus:
        _dist = lda_model.get_document_topics(_bow, minimum_probability=0)
        _distributions.append([_p for _, _p in sorted(_dist)])
        _dominant.append(max(_dist, key=lambda t: t[1])[0])

    topic_df = pd.DataFrame(
        _distributions,
        columns=[f"Topic_{_i}" for _i in range(optimal_num_topics)],
        index=otu_final.index,
    )
    dominant_topics = _dominant

    topic_df.head()
    return dominant_topics, topic_df


@app.cell
def _(dominant_topics, meta_final, topic_df):
    # Merge topics with metadata -
    results_df = meta_final.copy()
    results_df["Dominant_Topic"] = [f"Topic_{_t}" for _t in dominant_topics]

    for _col in topic_df.columns:
        results_df[_col] = topic_df[_col].values

    dominant_topic_counts = results_df["Dominant_Topic"].value_counts()
    print("Dominant topic distribution:")
    for _topic, _count in dominant_topic_counts.items():
        print(f"   {_topic}: {_count} samples")

    results_df.head()
    return dominant_topic_counts, results_df


@app.cell
def _(
    dictionary,
    lda_model,
    optimal_num_topics,
    otu_name_mapping,
    pd,
    taxonomy_filtered,
):
    # Taxa composition per topic 
    def get_top_taxa_for_topic(topic_id, top_n=10):
        """Return a DataFrame of the top taxa (with lineage) for one topic."""
        _rows = []
        for _term_id, _prob in lda_model.get_topic_terms(topic_id, topn=top_n):
            _term = dictionary[_term_id]
            _matches = [k for k, v in otu_name_mapping.items() if v == _term]
            if not _matches:
                continue
            _otu = _matches[0]
            if _otu not in taxonomy_filtered.index:
                continue
            _tax = taxonomy_filtered.loc[_otu]
            _rows.append(
                {
                    "term": _term,
                    "probability": _prob,
                    "phylum": _tax.get("Phylum", "Unknown"),
                    "class": _tax.get("Class", "Unknown"),
                    "order": _tax.get("Order", "Unknown"),
                    "family": _tax.get("Family", "Unknown"),
                    "genus": _tax.get("Genus", "Unknown"),
                }
            )
        return pd.DataFrame(_rows)

    _summaries = []
    for _topic_id in range(optimal_num_topics):
        _taxa = get_top_taxa_for_topic(_topic_id, top_n=20)
        if _taxa.empty:
            continue

        _phy = _taxa.groupby("phylum")["probability"].sum().sort_values(ascending=False)
        _fam = _taxa.groupby("family")["probability"].sum().sort_values(ascending=False)
        _gen = _taxa.groupby("genus")["probability"].sum().sort_values(ascending=False)

        _summaries.append(
            {
                "Topic_ID": f"Topic_{_topic_id}",
                "Top_Phylum_1": _phy.index[0] if len(_phy) else "Unknown",
                "Top_Phylum_1_Prob": _phy.values[0] if len(_phy) else 0,
                "Top_Phylum_2": _phy.index[1] if len(_phy) > 1 else "",
                "Top_Phylum_2_Prob": _phy.values[1] if len(_phy) > 1 else 0,
                "Top_Family_1": _fam.index[0] if len(_fam) else "Unknown",
                "Top_Family_1_Prob": _fam.values[0] if len(_fam) else 0,
                "Top_Genus_1": _gen.index[0] if len(_gen) else "Unknown",
                "Top_Genus_1_Prob": _gen.values[0] if len(_gen) else 0,
            }
        )

    topic_taxa_summary = pd.DataFrame(_summaries)
    topic_taxa_summary
    return (topic_taxa_summary,)


@app.cell
def _(dictionary, lda_model, optimal_num_topics, plt, sns, topic_taxa_summary):
    #Top taxa per topic
    _n_rows = (optimal_num_topics + 2) // 3
    _fig, _axes = plt.subplots(_n_rows, 3, figsize=(18, 4 * _n_rows))
    _axes = _axes.flatten()

    for _topic_id in range(optimal_num_topics):
        _ax = _axes[_topic_id]
        _terms = lda_model.get_topic_terms(_topic_id, topn=10)
        _words = [dictionary[_tid] for _tid, _ in _terms]
        _probs = [_p for _, _p in _terms]

        _row = topic_taxa_summary.loc[
            topic_taxa_summary["Topic_ID"] == f"Topic_{_topic_id}"
        ].iloc[0]
        _label = str(_row["Top_Phylum_1"])[:20]

        _ax.barh(_words, _probs, color=sns.color_palette("viridis", len(_words)))
        _ax.set_xlabel("Probability", fontsize=10)
        _display_id = _topic_id + 1
        _ax.set_title(
        f"Topic {_display_id} (MC{_display_id}) - {_label}",
        fontsize=12, fontweight="bold",)
        _ax.invert_yaxis()
        _ax.tick_params(axis="y", labelsize=8)

    for _idx in range(optimal_num_topics, len(_axes)):
        _axes[_idx].axis("off")

    _fig.tight_layout()
    _fig.savefig("topic_composition_with_taxa.png", dpi=300, bbox_inches="tight")
    _fig
    return


@app.cell
def _(
    dendrogram,
    linkage,
    mcolors,
    mpatches,
    np,
    optimal_num_topics,
    plt,
    results_df,
    sns,
    topic_df,
    topic_taxa_summary,
):
    #Clustered heatmap with metadata bars 
    def _make_palette(values, base_palette="tab20"):
        _vals = list(dict.fromkeys(values))
        _colors = sns.color_palette(base_palette, len(_vals))
        return {_v: mcolors.to_hex(_c) for _v, _c in zip(_vals, _colors)}

    def _draw_meta_bar(ax, labels, palette, title):
        _hex = [palette.get(_l, "#d3d3d3") for _l in labels]
        _rgba = np.array([mcolors.to_rgba(_c) for _c in _hex])[np.newaxis, :, :]
        ax.imshow(_rgba, aspect="auto")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_ylabel(title, rotation=0, ha="right", va="center", fontsize=10)

        _patches = [
            mpatches.Patch(color=palette.get(_l, "#d3d3d3"), label=_l)
            for _l in dict.fromkeys(labels)
        ]
        if _patches:
            ax.legend(
                handles=_patches,
                bbox_to_anchor=(1.01, 0.5),
                loc="center left",
                frameon=False,
                fontsize=8,
            )

    _col_linkage = linkage(topic_df.T.values, method="ward")
    _col_order = dendrogram(_col_linkage, no_plot=True)["leaves"]
    _row_linkage = linkage(topic_df.values, method="ward")
    _row_order = dendrogram(_row_linkage, no_plot=True)["leaves"]

    _ordered_labels = [f"Topic_{_i}" for _i in _col_order]
    _heat = topic_df.iloc[_row_order, _col_order].T.values

    _meta_cols = [
        _c for _c in ("Type", "Category", "Plastic") if _c in results_df.columns
    ]
    _n_meta = len(_meta_cols)

    _fig = plt.figure(figsize=(20, 12))
    _gs = _fig.add_gridspec(
        nrows=_n_meta + 2,
        ncols=2,
        height_ratios=[0.3] * _n_meta + [4, 0.3],
        width_ratios=[1, 0.05],
    )

    _meta_axes = [_fig.add_subplot(_gs[_i, 0]) for _i in range(_n_meta)]
    _ax_heat = _fig.add_subplot(_gs[_n_meta, 0])
    _ax_cbar = _fig.add_subplot(_gs[_n_meta, 1])
    _ax_dendro = _fig.add_subplot(_gs[_n_meta + 1, 0])

    _sample_order = topic_df.iloc[_row_order].index.tolist()
    _indexed = results_df.set_index("sample_name")

    for _i, _col in enumerate(_meta_cols):
        _vals = _indexed.loc[_sample_order, _col]
        _palette = _make_palette(_vals, base_palette="tab10" if _i == 0 else "Set2")
        _draw_meta_bar(_meta_axes[_i], _vals, _palette, _col)

    _im = _ax_heat.imshow(
        _heat, aspect="auto", cmap="Reds", interpolation="nearest", vmin=0, vmax=1
    )

    _y_labels = []
    for _label in _ordered_labels:
        _tid = int(_label.split("_")[1])

        # Lookup remains zero-based
        _row = topic_taxa_summary.loc[
            topic_taxa_summary["Topic_ID"] == f"Topic_{_tid}"
        ].iloc[0]

        # Display becomes one-based
        _display_id = _tid + 1
        _y_labels.append(
            f"MC{_display_id} ({str(_row['Top_Phylum_1'])[:15]})"
        )

    _ax_heat.set_yticks(range(len(_y_labels)))
    _ax_heat.set_yticklabels(_y_labels, fontsize=9)
    _ax_heat.set_xticks([])
    _ax_heat.set_ylabel("Microbiome Components (dominant taxa)", fontsize=11)

    _fig.colorbar(_im, cax=_ax_cbar, label="MC Probability")

    dendrogram(
        _col_linkage,
        ax=_ax_dendro,
        labels=[f"MC{_i + 1}" for _i in range(optimal_num_topics)],
        orientation="bottom",
        leaf_font_size=8,
    )
    _ax_dendro.set_xlabel("Distance")
    _ax_dendro.set_ylabel("")

    _fig.suptitle(
        "Microbiome Topic Modeling with Taxa Annotations",
        fontsize=16, fontweight="bold",
    )
    _fig.tight_layout()
    _fig.savefig("lda_heatmap_with_taxa_labels.png", dpi=300, bbox_inches="tight")
    _fig
    return


@app.cell
def _(optimal_num_topics, plt, results_df, sns):
    # topic probability by sample Type
    if "Type" in results_df.columns:
        _n_rows = (optimal_num_topics + 2) // 3
        _fig, _axes = plt.subplots(_n_rows, 3, figsize=(18, 4 * _n_rows))
        _axes = _axes.flatten()

        for _topic_id in range(optimal_num_topics):
            _ax = _axes[_topic_id]
            sns.boxplot(
                data=results_df, x="Type", y=f"Topic_{_topic_id}",
                ax=_ax, palette="Set2",
            )
            _display_id = _topic_id + 1
            _ax.set_title(
            f"Topic {_display_id} by Type",
                fontsize=11,
                fontweight="bold")
            _ax.set_ylabel(
                f"MC{_display_id} Probability",
                fontsize=10
            )
        
            _ax.set_ylabel(f"MC{_topic_id} Probability", fontsize=10)
            _ax.set_xlabel("Type", fontsize=10)
            _ax.tick_params(axis="x", rotation=45)

        for _idx in range(optimal_num_topics, len(_axes)):
            _axes[_idx].axis("off")

        _fig.tight_layout()
        _fig.savefig("topic_distributions_by_type.png", dpi=300, bbox_inches="tight")
        _out = _fig
    else:
        _out = "No 'Type' column in metadata - box plots skipped."

    _out
    return


@app.cell
def _(
    TOPIC_MIN,
    best_coherence,
    coherence_values,
    dictionary,
    dominant_topic_counts,
    lda_model,
    optimal_num_topics,
    otu_final,
    perplexity_values,
    results_df,
    topic_taxa_summary,
):
    #Persist model, tables and a text report 
    lda_model.save("lda_model")
    dictionary.save("lda_dictionary")
    results_df.to_csv("final_results_with_metadata.csv", index=False)
    topic_taxa_summary.to_csv("topic_taxa_summary.csv", index=False)

    with open("model_statistics.txt", "w") as _f:
        _f.write("=" * 80 + "\n")
        _f.write("LDA MODEL STATISTICS WITH TAXA INFORMATION\n")
        _f.write("=" * 80 + "\n\n")
        _f.write(f"Samples (after filtering): {otu_final.shape[0]}\n")
        _f.write(f"OTUs (after filtering):    {otu_final.shape[1]}\n")
        _f.write(f"Number of topics:          {optimal_num_topics}\n")
        _f.write(f"Coherence score:           {best_coherence:.4f}\n")
        _f.write(
            f"Log perplexity:            "
            f"{perplexity_values[optimal_num_topics - TOPIC_MIN]:.4f}\n"
        )

        _f.write("\n" + "=" * 80 + "\nTAXA COMPOSITION OF EACH TOPIC\n" + "=" * 80 + "\n")
        for _, _row in topic_taxa_summary.iterrows():
            _f.write(f"\n{_row['Topic_ID']}:\n")
            _f.write(
                f"  Primary Phylum:   {_row['Top_Phylum_1']} "
                f"({_row['Top_Phylum_1_Prob']:.3f})\n"
            )
            if _row["Top_Phylum_2"]:
                _f.write(
                    f"  Secondary Phylum: {_row['Top_Phylum_2']} "
                    f"({_row['Top_Phylum_2_Prob']:.3f})\n"
                )
            _f.write(
                f"  Primary Family:   {_row['Top_Family_1']} "
                f"({_row['Top_Family_1_Prob']:.3f})\n"
            )
            _f.write(
                f"  Primary Genus:    {_row['Top_Genus_1']} "
                f"({_row['Top_Genus_1_Prob']:.3f})\n"
            )

        _f.write("\n" + "=" * 80 + "\nSAMPLE DISTRIBUTION\n" + "=" * 80 + "\n")
        for _col in ("Category", "Type", "Plastic"):
            if _col in results_df.columns:
                _f.write(f"\nBy {_col}:\n")
                for _val, _count in results_df[_col].value_counts().items():
                    _f.write(f"  {_val}: {_count}\n")

        _f.write("\n" + "=" * 80 + "\nDOMINANT TOPIC DISTRIBUTION\n" + "=" * 80 + "\n")
        for _topic, _count in dominant_topic_counts.items():
            _f.write(f"  {_topic}: {_count} samples\n")

    print("Saved:")
    for _name in (
        "model_selection_metrics.png",
        "topic_composition_with_taxa.png",
        "lda_heatmap_with_taxa_labels.png",
        "topic_distributions_by_type.png",
        "final_results_with_metadata.csv",
        "topic_taxa_summary.csv",
        "model_statistics.txt",
        "lda_model",
        "lda_dictionary",
    ):
        print(f"  - {_name}")
    _ = coherence_values 
    return


@app.cell
def _():
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
