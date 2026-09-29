import marimo

__generated_with = "0.23.10"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import numpy as np

    from svg import SVG, Rect, G, Text, Title

    from scipy.cluster.hierarchy import linkage, dendrogram

    import warnings
    warnings.filterwarnings("ignore")
    return G, Rect, SVG, Text, Title, dendrogram, linkage, mo, np, pd


@app.cell
def _(pd):
    df = pd.read_csv("final_results_with_metadata.csv")
    if "sample_name" in df.columns:
        df = df.set_index("sample_name")



    # rename exposure time names
    CATEGORY_RENAME = {
        "Daphnia: Pre-exposure : Clone control": "Daphnia: Pre-exposure-Adult",
        "Daphnia: Juvenile: Post-exposure": "Daphnia: 23-Day exposure-Juvenile",
        "Daphnia: Post-exposure": "Daphnia: 23-Day exposure-Adult",
        "Daphnia: Post-Long-exposure (selected survivals)": "Daphnia: Long-exposure survivors",
    }
    if "Category" in df.columns:
        _norm = df["Category"].astype(str).str.strip().str.replace(r"\s+", " ", regex=True)
        df["Category"] = _norm.replace(CATEGORY_RENAME)

    topic_cols = [c for c in df.columns if c.startswith("Topic_")]
    topic_df = df[topic_cols].copy()
    print(f"LDA topic matrix: {topic_df.shape[0]} samples x {topic_df.shape[1]} topics")
    return df, topic_df


@app.cell
def _(df):
    df
    return


@app.cell
def _(df):
    TRACK_ORDER = ["Pond", "Plastic", "Category","Clone"]

    meta_cols = [c for c in TRACK_ORDER if c in df.columns]
    meta_df = df[meta_cols].copy()

    # fillna 
    for c in meta_cols:
        meta_df[c] = meta_df[c].fillna("n/a").astype(str).str.strip()
        meta_df[c] = meta_df[c].replace(
            {"nan": "n/a", "NaN": "n/a", "None": "n/a", "": "n/a"}
        )

    meta_df = meta_df.reset_index().rename(columns={"index": "sample_name"})
    print("Annotation tracks used:", meta_cols)
    return TRACK_ORDER, meta_cols, meta_df


@app.cell
def _(meta_df):
    meta_df
    return


@app.cell
def _(np, topic_df):
    sample_counts = topic_df.sum(axis=1)
    count_25th = np.percentile(sample_counts, 25)
    count_50th = np.percentile(sample_counts, 50)
    count_75th = np.percentile(sample_counts, 75)
    return count_25th, count_50th, count_75th, sample_counts


@app.cell
def _(dendrogram, linkage, meta_cols, meta_df, np, topic_df):
    # Sort orders.
    # "By Dominant Topic": group samples by strongest topic 
    #"Hierarchical Clustering": Ward linkage on the full topic profile 
    #groups samples by overall similarity across all topics.
    sort_orders = {}

    _tv = topic_df.values
    sort_orders["By Dominant Topic"] = list(np.lexsort((-_tv.max(axis=1), _tv.argmax(axis=1))))

    if topic_df.shape[0] > 1:
        _link = linkage(topic_df.values, method="ward")
        sort_orders["Hierarchical Clustering"] = dendrogram(_link, no_plot=True)["leaves"]

    sort_orders["Original Order"] = list(range(len(meta_df)))

    _df = meta_df.copy()
    _df["idx"] = np.arange(len(_df))
    for _col in meta_cols:
        sort_orders[f"By {_col}"] = _df.sort_values([_col, "idx"]).index.tolist()

    print("Sort orders:", list(sort_orders.keys()))
    return (sort_orders,)


@app.cell
def _(meta_df, mo, sort_orders):
    options = list(sort_orders.keys())
    sort_dropdown = mo.ui.dropdown(
        options=options, value="By Dominant Topic", label="Sort samples by:"
    )
    count_filter = mo.ui.dropdown(
        options=["All samples", "Top 75%", "Top 50%", "Top 25%"],
        value="All samples", label="Filter by total count:",
    )

    filter_specs = [
        ("type_filter", "Type", "Type:"),
        ("pond_filter", "Pond", "Inoculum Pond:"),
        ("clone_filter", "Clone", "Host clone:"),
        ("plastic_filter", "Plastic", "Plastic Types:"),
        ("category_filter", "Category", "Exposure Category:"),
    ]
    filter_dropdowns = {}
    for _key, _col, _lbl in filter_specs:
        if _col in meta_df.columns:
            _opts = ["All"] + sorted(meta_df[_col].unique().tolist())
            filter_dropdowns[_key] = mo.ui.dropdown(options=_opts, value="All", label=_lbl)

    filter_col = {
        "type_filter": "Type", "pond_filter": "Pond",
        "clone_filter": "Clone", 
        "plastic_filter": "Plastic", "category_filter": "Category",
    }

    mo.hstack([sort_dropdown, count_filter] + list(filter_dropdowns.values()), justify="start")
    return count_filter, filter_col, filter_dropdowns, sort_dropdown


@app.cell
def _(
    G,
    Rect,
    SVG,
    TRACK_ORDER,
    Text,
    Title,
    dendrogram,
    filter_col,
    linkage,
    meta_df,
    np,
    topic_df,
):
    FONT = "Helvetica, Arial, sans-serif"

    GREY_LIGHT, GREY_MID, GREY_DARK, GREY_NA = "#e0e0e0", "#9e9e9e", "#616161", "#f2f2f2"
    BLUE, SKY, GREEN, PURPLE, ORANGE, YELLOW, TEAL = (
        "#0072B2", "#56B4E9", "#009E73", "#CC79A7", "#E69F00", "#F0E442", "#1B9E9E"
    )

    TRACK_PALETTES = {
        "Type": {"Bacterioplankton": BLUE, "Daphnia": GREEN},
        "Pond": {"BP": SKY, "DG": PURPLE, "H2O": GREY_MID},
        "Clone": {
            "BH": BLUE, "F": GREEN,
            "KNO15.04": ORANGE, "KNO15": "#B36B00", "KNO15.05": "#F4C771",
            "C": GREY_LIGHT, "Not sure": GREY_MID, "H2O": GREY_DARK,
        },
        "Plastic": {
            "Nylon": BLUE, "PET": GREEN, "PLA": PURPLE,
            "C": GREY_LIGHT, "Control": GREY_MID, "Not sure": GREY_DARK,
        },
        "Category": {
            "Daphnia: Pre-exposure-Adult": YELLOW,
            "Daphnia: 23-Day exposure-Juvenile": ORANGE,
            "Daphnia: 23-Day exposure-Adult": BLUE,
            "Daphnia: Long-exposure survivors": PURPLE,
        },
    }
    TRACK_TITLE = {
        "Type": "Type",
        "Pond": "Inoculum Pond",
        "Clone": "Host Clone",
        "Plastic": "Plastic Types",
        "Category": "Exposure Category",
    }
    _FALLBACK = [TEAL, "#984ea3", "#8c8c8c", "#bcbd22", "#17becf"]

    def meta_to_color(label, track):
        if label in ("n/a", "nan", "None", ""):
            return GREY_NA
        pal = TRACK_PALETTES.get(track, {})
        return pal.get(label) or _FALLBACK[abs(hash(label)) % len(_FALLBACK)]

    PROB_COLORS = ["#fff5f0", "#fee0d2", "#fcbba1", "#fc9272",
                   "#fb6a4a", "#ef3b2c", "#cb181d", "#99000d"]

    def prob_to_color(p):
        if p is None or (isinstance(p, float) and np.isnan(p)):
            return "#ffffff"
        p = max(0.0, min(1.0, float(p)))
        return PROB_COLORS[int(round(p * (len(PROB_COLORS) - 1)))]

    def T(**kw):
        kw.setdefault("font_family", FONT)
        return Text(**kw)


    #create a heatmap
    def create_lda_heatmap(
        sort_method="By Dominant Topic", sort_orders=None,
        count_filter_value="All samples", filter_values=None,
        sample_counts=None, count_thresholds=None,
    ):
        if filter_values is None:
            filter_values = {}

        samples = meta_df.copy().reset_index(drop=True)
        topics = topic_df.copy().reset_index(drop=True)
        tcols = [c for c in topics.columns if c.startswith("Topic_")]

        #filtering
        mask = np.ones(len(samples), dtype=bool)
        if count_filter_value != "All samples" and sample_counts is not None:
            thr = {"Top 75%": count_thresholds[0], "Top 50%": count_thresholds[1],
                   "Top 25%": count_thresholds[2]}[count_filter_value]
            mask &= (sample_counts.values >= thr)
        for fkey, col in filter_col.items():
            val = filter_values.get(fkey, "All")
            if val != "All" and col in samples.columns:
                mask &= (samples[col] == val).values
        samples = samples[mask].reset_index(drop=True)
        topics = topics[mask].reset_index(drop=True)

        if len(samples) == 0:
            return SVG(width=520, height=90, elements=[
                T(x=260, y=50, text="No samples match the selected filters",
                  text_anchor="middle", font_size="14", fill="#cc0000")])

        # sorting
        tv = topics[tcols].values
        if sort_method == "By Dominant Topic":
            order = list(np.lexsort((-tv.max(axis=1), tv.argmax(axis=1))))
        elif sort_method == "Hierarchical Clustering" and len(samples) > 1:
            order = dendrogram(linkage(tv, method="ward"), no_plot=True)["leaves"]
        elif sort_method.startswith("By ") and sort_method[3:] in samples.columns:
            col = sort_method[3:]
            samples["_i"] = np.arange(len(samples))
            order = samples.sort_values([col, "_i"]).index.tolist()
            samples = samples.drop("_i", axis=1)
        else:
            order = list(range(len(samples)))
        samples = samples.iloc[order].reset_index(drop=True)
        topics = topics.iloc[order].reset_index(drop=True)

        #geometry 
        n = len(samples)
        n_topics = len(tcols)
        tracks = [t for t in TRACK_ORDER if t in samples.columns]

        col_w = 5.5                    # column pitch per sample (less compact)
        cell_w = 3.4                   # bar narrower than pitch -> visible gap
        cell_off = (col_w - cell_w) / 2  # centre each bar in its slot
        ann_h, ann_vgap = 16.0, 3.0
        mc_h, mc_vgap = 26.0, 2.5
        gap_ann_heat = 18.0
        y_title_x = 14
        label_w = 86                   # right-aligned track / MC labels
        x0 = y_title_x + label_w
        right_pad = 28
        plot_w = n * col_w
        width = x0 + plot_w + right_pad

        elements = []
        elements.append(T(x=x0, y=18, text="Microbial Community LDA",
                          font_size="14", font_weight="bold", fill="#222"))
        elements.append(T(x=x0, y=31, fill="#888", font_size="9",
                          text=f"{n} samples \u00b7 {n_topics} topics \u00b7 sort: {sort_method}"))

        y = 44
        # annotation tracks 
        for track in tracks:
            elements.append(T(x=x0 - 6, y=y + ann_h - 4, text=TRACK_TITLE.get(track, track),
                              text_anchor="end", font_size="9", fill="#333"))
            vals = samples[track].tolist()
            for j, val in enumerate(vals):
                elements.append(Rect(
                    x=x0 + j * col_w + cell_off, y=y, width=cell_w, height=ann_h,
                    fill=meta_to_color(val, track), stroke="none",
                    elements=[Title(elements=[f"{TRACK_TITLE.get(track, track)}: {val}"])]))
            y += ann_h + ann_vgap

        y += gap_ann_heat - ann_vgap
        heat_top = y

        # topic heatmap 
        for i, tcol in enumerate(tcols):
            elements.append(T(x=x0 - 6, y=y + mc_h / 2 + 3, text=f"MC{i + 1}",
                              text_anchor="end", font_size="9.5", fill="#333"))
            sn = samples["sample_name"].tolist()
            pv = topics[tcol].tolist()
            for j in range(n):
                elements.append(Rect(
                    x=x0 + j * col_w + cell_off, y=y, width=cell_w, height=mc_h,
                    fill=prob_to_color(pv[j]), stroke="none",
                    elements=[Title(elements=[f"{sn[j]} \u00b7 MC{i + 1}: p={float(pv[j]):.3f}"])]))
            y += mc_h + mc_vgap
        heat_bottom = y - mc_vgap

        # axis titles 
        elements.append(G(
            elements=[T(x=0, y=0, text="Microbial Components", text_anchor="middle",
                        font_size="11", fill="#333")],
            transform=f"translate({y_title_x},{(heat_top + heat_bottom)/2}) rotate(-90)"))
        elements.append(T(x=x0 + plot_w / 2, y=heat_bottom + 20, text="Samples",
                          text_anchor="middle", font_size="11", fill="#333"))


        #  LEGENDS  

        sw, lh, cw = 11, 15, 6.0      # swatch, line height, char width estimate
        gx0 = x0
        y = heat_bottom + 46
        x = gx0
        row_max_h = 0
        gutter = 34
        max_x = width - 10

        def group_size(title, items):
            label_lens = [len(str(t)) for t in items]
            long = max(label_lens) if label_lens else 0
            ncol = 1 if long > 20 else 2
            rows = (len(items) + ncol - 1) // ncol
            col_w_px = sw + 6 + long * cw + 18
            gw = max(len(title) * 6.5, ncol * col_w_px)
            gh = 16 + rows * lh
            return gw, gh, ncol, col_w_px

        def draw_group(gx, gy, title, items, track):
            gw, gh, ncol, col_w_px = group_size(title, items)
            elements.append(T(x=gx, y=gy + 10, text=title, font_size="10",
                              font_weight="bold", fill="#222"))
            yy = gy + 16
            for k, val in enumerate(items):
                ci = k % ncol
                ri = k // ncol
                ix = gx + ci * col_w_px
                iy = yy + ri * lh
                elements.append(Rect(x=ix, y=iy, width=sw, height=sw,
                                     fill=meta_to_color(val, track),
                                     stroke="#888", stroke_width="0.5"))
                elements.append(T(x=ix + sw + 5, y=iy + sw - 1, text=str(val),
                                  font_size="8.5", fill="#444"))
            return gw, gh

        # one legend group per categorical track
        for track in tracks:
            present = list(dict.fromkeys(samples[track].tolist()))
            pal = TRACK_PALETTES.get(track, {})
            items = [v for v in pal if v in present] + \
                    [v for v in present if v not in pal and v != "n/a"]
            if "n/a" in present:
                items.append("n/a")
            gw, gh, _, _ = group_size(TRACK_TITLE.get(track, track), items)
            if x + gw > max_x and x > gx0:
                x = gx0
                y += row_max_h + 14
                row_max_h = 0
            _, dh = draw_group(x, y, TRACK_TITLE.get(track, track), items, track)
            x += gw + gutter
            row_max_h = max(row_max_h, dh)

        #MC probability gradient group 
        grad_w = 150
        gtitle, gh = "MC Probability", 16 + lh + 12
        if x + grad_w > max_x and x > gx0:
            x = gx0
            y += row_max_h + 14
            row_max_h = 0
        elements.append(T(x=x, y=y + 10, text=gtitle, font_size="10",
                          font_weight="bold", fill="#222"))
        seg = grad_w / len(PROB_COLORS)
        gy = y + 16
        for i, c in enumerate(PROB_COLORS):
            elements.append(Rect(x=x + i * seg, y=gy, width=seg + 0.6, height=sw,
                                 fill=c, stroke="none"))
        elements.append(Rect(x=x, y=gy, width=grad_w, height=sw,
                             fill="none", stroke="#888", stroke_width="0.5"))
        for i, v in enumerate([0, 0.5, 1.0]):
            elements.append(T(x=x + i * (grad_w / 2), y=gy + sw + 9, text=f"{v:.1f}",
                              font_size="7.5", fill="#777", text_anchor="middle"))
        row_max_h = max(row_max_h, gh)

        height = y + row_max_h + 16
        return SVG(width=width, height=height, elements=elements)

    return (create_lda_heatmap,)


@app.cell
def _(
    count_25th,
    count_50th,
    count_75th,
    count_filter,
    create_lda_heatmap,
    filter_dropdowns,
    mo,
    sample_counts,
    sort_dropdown,
    sort_orders,
):
    svg = create_lda_heatmap(
        sort_dropdown.value, sort_orders=sort_orders,
        count_filter_value=count_filter.value,
        filter_values={k: d.value for k, d in filter_dropdowns.items()},
        sample_counts=sample_counts,
        count_thresholds=(count_25th, count_50th, count_75th),
    )
    ui = [sort_dropdown, count_filter] + list(filter_dropdowns.values())
    mo.vstack([mo.hstack(ui, justify="start"), mo.Html(svg.as_str())])
    return


if __name__ == "__main__":
    app.run()
