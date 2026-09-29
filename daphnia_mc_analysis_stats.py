import marimo

__generated_with = "0.23.10"
app = marimo.App(width="medium")


@app.cell
def _():
    import re
    import warnings

    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from scipy import stats
    from scipy.spatial.distance import jensenshannon

    warnings.filterwarnings("ignore")
    return jensenshannon, mo, np, pd, plt, re, stats


@app.cell
def _(pd, re):

    # Load and prepare data

    df = pd.read_csv("final_results_with_metadata.csv")

    if "sample_name" in df.columns:
        df = df.set_index("sample_name")

    # Daphnia samples only
    df = df[
        df["Type"].astype(str).str.strip().str.casefold() == "daphnia"
    ].copy()

    # Keep only columns named Topic_0, Topic_1, ...
    topic_cols = [
        c for c in df.columns
        if re.fullmatch(r"Topic_\d+", str(c))
    ]

    # Sort numerically, so Topic_10 does not appear before Topic_2
    topic_cols = sorted(
        topic_cols,
        key=lambda c: int(re.search(r"\d+$", c).group())
    )

    if not topic_cols:
        raise ValueError(
            "No topic columns were found. Expected names such as "
            "'Topic_0', 'Topic_1', ..."
        )

    topic_labels = [
        f"MC{i}"
        for i in range(len(topic_cols))
    ]

    P = df[topic_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)

    # Remove samples whose topic probabilities do not have a positive sum
    valid_rows = P.sum(axis=1) > 0
    if not valid_rows.all():
        print(
            f"Removed {(~valid_rows).sum()} sample(s) with zero or invalid "
            "total topic probability."
        )
        df = df.loc[valid_rows].copy()
        P = P.loc[valid_rows].copy()

    # Renormalise each sample so topic proportions sum exactly to one
    P = P.div(P.sum(axis=1), axis=0)

    candidate_vars = [
        "Pond",
        "Clone",
        "Plastic",
        "Category",
        "Plastic_content",
        "Comparison",
    ]

    meta = df[
        [v for v in candidate_vars if v in df.columns]
    ].copy()

    for c in meta.columns:
        meta[c] = (
            meta[c]
            .fillna("n/a")
            .astype(str)
            .str.strip()
            .replace({
                "nan": "n/a",
                "NaN": "n/a",
                "None": "n/a",
                "": "n/a",
            })
        )

    def _testable(col):
        _vc = meta.loc[meta[col] != "n/a", col].value_counts()
        return int((_vc >= 3).sum()) >= 2

    meta_vars = [c for c in meta.columns if _testable(c)]
    meta = meta[meta_vars]

    print(
        f"DAPHNIA ONLY · {P.shape[0]} samples · "
        f"{len(topic_cols)} microbial components"
    )
    print("Topic columns:", topic_cols)
    print("Metadata variables tested:", meta_vars)


    return P, meta, meta_vars, topic_cols, topic_labels


@app.cell
def _(np):

    # Statistical analysis

    def bh_fdr(pvals):
        """Benjamini-Hochberg FDR q-values aligned to input order."""
        p = np.asarray(pvals, dtype=float)

        if p.size == 0:
            return np.asarray([], dtype=float)

        n = p.size
        order = np.argsort(p)
        ranked = p[order] * n / (np.arange(n) + 1)

        # Enforce monotonicity from the largest p-value downward
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]

        q = np.empty(n, dtype=float)
        q[order] = np.clip(ranked, 0.0, 1.0)
        return q

    def kruskal_eta2(H, k, n):
        """
        Bias-adjusted Kruskal-Wallis effect-size estimate:
        (H - k + 1) / (n - k).

        Small negative values can occur under the null hypothesis.
        They are retained in the results table but clipped to zero
        only when plotting.
        """
        denominator = n - k
        if denominator <= 0:
            return np.nan
        return float((H - k + 1) / denominator)

    def clr(mat, pseudo=None):
        """Centred log-ratio transform with simple zero replacement."""
        X = np.asarray(mat, dtype=float)

        if pseudo is None:
            nonzero = X[X > 0]
            pseudo = nonzero.min() / 2.0 if nonzero.size else 1e-6

        X = np.where(X <= 0, pseudo, X)
        X = X / X.sum(axis=1, keepdims=True)

        logX = np.log(X)
        return logX - logX.mean(axis=1, keepdims=True)

    def significance_stars(p):
        """plot legend."""
        if p <= 0.001:
            return "***"
        if p <= 0.01:
            return "**"
        if p <= 0.05:
            return "*"
        return ""

    return bh_fdr, clr, kruskal_eta2, significance_stars


@app.cell
def _(P, bh_fdr, kruskal_eta2, meta, meta_vars, pd, stats, topic_cols):

    # Component-wise Kruskal-Wallis screening

    def _screen(min_group=3):
        _rows = []

        for _var in meta_vars:
            _labels = meta[_var]

            for _tcol in topic_cols:
                _sub = pd.DataFrame({
                    "y": P[_tcol],
                    "g": _labels,
                }).dropna()

                _sub = _sub[_sub["g"] != "n/a"]

                _counts = _sub["g"].value_counts()
                _keep = _counts[_counts >= min_group].index
                _sub = _sub[_sub["g"].isin(_keep)]

                _k = _sub["g"].nunique()
                if _k < 2:
                    continue

                _groups = [
                    _group["y"].to_numpy()
                    for _, _group in _sub.groupby("g", sort=True)
                ]

                _H, _p = stats.kruskal(*_groups)

                _rows.append({
                    "variable": _var,
                    "topic": _tcol,
                    "n_groups": int(_k),
                    "n": int(len(_sub)),
                    "H": float(_H),
                    "p": float(_p),
                    "eta2": kruskal_eta2(_H, _k, len(_sub)),
                })

        _out = pd.DataFrame(_rows)

        if _out.empty:
            return _out

        # Global BH correction across the complete variable × component grid
        _out["q_fdr"] = bh_fdr(_out["p"].to_numpy())

        return _out.sort_values(
            ["q_fdr", "variable", "topic"]
        ).reset_index(drop=True)

    results = _screen()

    if results.empty:
        print("No valid component-wise tests were available.")
    else:
        print(
            f"{len(results)} tests · "
            f"{(results['q_fdr'] < 0.05).sum()} significant at global FDR < 0.05"
        )
    return (results,)


@app.cell
def _(results):
    results.round(4)
    return


@app.cell
def _(
    meta_vars,
    np,
    plt,
    results,
    significance_stars,
    topic_cols,
    topic_labels,
):

    # Effect-size heatmap

    def _heatmap():
        _eta_raw = np.full(
            (len(topic_cols), len(meta_vars)),
            np.nan,
            dtype=float,
        )
        _ann = np.full(
            (len(topic_cols), len(meta_vars)),
            "",
            dtype=object,
        )

        _topic_index = {topic: i for i, topic in enumerate(topic_cols)}
        _variable_index = {
            variable: j for j, variable in enumerate(meta_vars)
        }

        for _, _row in results.iterrows():
            _i = _topic_index[_row["topic"]]
            _j = _variable_index[_row["variable"]]

            _eta_raw[_i, _j] = _row["eta2"]
            _ann[_i, _j] = significance_stars(_row["q_fdr"])

        # Negative estimates mean approximately zero effect.
        # Clip them only for visualisation.
        _eta_plot = np.clip(_eta_raw, 0.0, None)

        _finite = _eta_plot[np.isfinite(_eta_plot)]
        _vmax = float(_finite.max()) if _finite.size else 1.0
        if _vmax <= 0:
            _vmax = 1.0

        _fig, _ax = plt.subplots(
            figsize=(
                1.25 * len(meta_vars) + 2.2,
                0.55 * len(topic_cols) + 2.2,
            )
        )

        _im = _ax.imshow(
            _eta_plot,
            cmap="PuRd",
            vmin=0,
            vmax=_vmax,
            aspect="auto",
        )

        _ax.set_xticks(range(len(meta_vars)))
        _ax.set_xticklabels(meta_vars, rotation=30, ha="right")

        _ax.set_yticks(range(len(topic_cols)))
        _ax.set_yticklabels(topic_labels)

        for _i in range(len(topic_cols)):
            for _j in range(len(meta_vars)):
                if np.isnan(_eta_raw[_i, _j]):
                    continue

                _display_value = max(0.0, _eta_raw[_i, _j])
                _star = _ann[_i, _j]
                _label = (
                    f"{_display_value:.2f}\n{_star}"
                    if _star
                    else f"{_display_value:.2f}"
                )

                _ax.text(
                    _j,
                    _i,
                    _label,
                    ha="center",
                    va="center",
                    fontsize=8,
                    color=(
                        "white"
                        if _display_value > _vmax * 0.55
                        else "#333333"
                    ),
                )

        _ax.set_title(
            "Effect size η² (Kruskal-Wallis)"
            "  ·  * q≤.05  ** q≤.01  *** q≤.001"
        )

        _fig.colorbar(
            _im,
            ax=_ax,
            fraction=0.04,
            pad=0.02,
            label="η²",
        )

        _fig.tight_layout()
        return _fig

    _heatmap()
    return


@app.cell
def _(np):

    # PERMANOVA and PERMDISP helpers

    def _ss_within(D2, codes):
        """Within-group sum of squares from a squared-distance matrix."""
        _ssw = 0.0

        for _group in np.unique(codes):
            _idx = np.where(codes == _group)[0]
            _n_group = len(_idx)

            if _n_group > 1:
                _block = D2[np.ix_(_idx, _idx)]
                _ssw += _block.sum() / (2.0 * _n_group)

        return float(_ssw)

    def permanova(D2, codes, perms=9999, seed=0):
        """
        One-factor PERMANOVA using unrestricted label permutations.

        D2 must be a squared-distance matrix.
        """
        _codes = np.asarray(codes)
        _n = len(_codes)
        _n_groups = len(np.unique(_codes))

        if _n_groups < 2:
            raise ValueError("PERMANOVA requires at least two groups.")

        _ss_total = D2.sum() / (2.0 * _n)
        _ss_within_observed = _ss_within(D2, _codes)
        _ss_among = _ss_total - _ss_within_observed

        _df_among = _n_groups - 1
        _df_within = _n - _n_groups

        _F_observed = (
            (_ss_among / _df_among)
            / (_ss_within_observed / _df_within)
        )
        _R2 = _ss_among / _ss_total

        _rng = np.random.default_rng(seed)
        _greater_equal = 1

        for _ in range(perms):
            _permuted_codes = _rng.permutation(_codes)
            _ss_within_perm = _ss_within(D2, _permuted_codes)
            _ss_among_perm = _ss_total - _ss_within_perm

            _F_perm = (
                (_ss_among_perm / _df_among)
                / (_ss_within_perm / _df_within)
            )

            if _F_perm >= _F_observed:
                _greater_equal += 1

        _p_value = _greater_equal / (perms + 1)

        return (
            float(_F_observed),
            float(_R2),
            float(_p_value),
        )

    def pcoa_coordinates(D2, tolerance=1e-10):
        """
        Principal coordinates from a squared-distance matrix.

        Positive-coordinate axes are returned. For the current topic-distance
        analyses this provides the Euclidean coordinate space used to calculate
        distances to group centroids.
        """
        _n = D2.shape[0]
        _J = np.eye(_n) - np.ones((_n, _n)) / _n
        _B = -0.5 * _J @ D2 @ _J

        _eigenvalues, _eigenvectors = np.linalg.eigh(_B)
        _order = np.argsort(_eigenvalues)[::-1]
        _eigenvalues = _eigenvalues[_order]
        _eigenvectors = _eigenvectors[:, _order]

        _positive = _eigenvalues > tolerance

        if not np.any(_positive):
            raise ValueError(
                "PCoA produced no positive eigenvalues for this distance matrix."
            )

        _coordinates = (
            _eigenvectors[:, _positive]
            * np.sqrt(_eigenvalues[_positive])
        )

        return _coordinates, _eigenvalues

    def _dispersion_f(coordinates, codes):
        """
        ANOVA F statistic for distances from observations to their
        group centroids in PCoA space.
        """
        _codes = np.asarray(codes)
        _distances = np.zeros(len(_codes), dtype=float)

        for _group in np.unique(_codes):
            _idx = np.where(_codes == _group)[0]
            _centroid = coordinates[_idx].mean(axis=0)

            _distances[_idx] = np.sqrt(
                ((coordinates[_idx] - _centroid) ** 2).sum(axis=1)
            )

        _grand_mean = _distances.mean()
        _ss_between = 0.0
        _ss_within = 0.0

        for _group in np.unique(_codes):
            _group_distances = _distances[_codes == _group]
            _group_mean = _group_distances.mean()

            _ss_between += len(_group_distances) * (
                _group_mean - _grand_mean
            ) ** 2

            _ss_within += np.sum(
                (_group_distances - _group_mean) ** 2
            )

        _n = len(_distances)
        _k = len(np.unique(_codes))

        _df_between = _k - 1
        _df_within = _n - _k

        _ms_between = _ss_between / _df_between
        _ms_within = _ss_within / _df_within

        if _ms_within <= 0:
            _F = np.inf if _ms_between > 0 else 0.0
        else:
            _F = _ms_between / _ms_within

        return float(_F), _distances

    def permdisp(D2, codes, perms=9999, seed=0):
        """
        Permutation test for equality of multivariate dispersion.

        Group labels are permuted and distances to the resulting group
        centroids are recomputed for every permutation.
        """
        _codes = np.asarray(codes)
        _coordinates, _ = pcoa_coordinates(D2)

        _F_observed, _distances = _dispersion_f(
            _coordinates,
            _codes,
        )

        _rng = np.random.default_rng(seed)
        _greater_equal = 1

        for _ in range(perms):
            _permuted_codes = _rng.permutation(_codes)
            _F_perm, _ = _dispersion_f(
                _coordinates,
                _permuted_codes,
            )

            if _F_perm >= _F_observed:
                _greater_equal += 1

        _p_value = _greater_equal / (perms + 1)

        return (
            float(_F_observed),
            float(_p_value),
            _distances,
        )

    return permanova, permdisp


@app.cell
def _(P, clr, jensenshannon, meta, meta_vars, np, pd, permanova, permdisp):

    # Run one-factor PERMANOVA and PERMDISP

    def _run_permanova(min_group=3, perms=9999):
        from scipy.spatial.distance import pdist as _pdist
        from scipy.spatial.distance import squareform as _squareform

        _P_normalised = P.div(P.sum(axis=1), axis=0).to_numpy()

        # Aitchison distance = Euclidean distance in CLR space
        _clr_values = clr(_P_normalised)
        _D_aitchison = _squareform(
            _pdist(_clr_values, metric="euclidean")
        )

        # Jensen-Shannon distance
        _n = len(_P_normalised)
        _D_js = np.zeros((_n, _n), dtype=float)

        for _i in range(_n):
            for _j in range(_i + 1, _n):
                _distance = jensenshannon(
                    _P_normalised[_i],
                    _P_normalised[_j],
                )
                _distance = (
                    0.0 if np.isnan(_distance) else float(_distance)
                )
                _D_js[_i, _j] = _distance
                _D_js[_j, _i] = _distance

        _rows = []

        for _var_index, _var in enumerate(meta_vars):
            _labels = meta[_var].to_numpy()

            _nonmissing = _labels != "n/a"
            _counts = pd.Series(
                _labels[_nonmissing]
            ).value_counts()

            _keep = _counts[_counts >= min_group].index

            _keep_mask = (
                _nonmissing
                & pd.Series(_labels).isin(_keep).to_numpy()
            )

            _sample_indices = np.where(_keep_mask)[0]
            _kept_labels = _labels[_keep_mask]
            _codes = pd.Categorical(_kept_labels).codes

            if len(np.unique(_codes)) < 2:
                continue

            for _distance_index, (_name, _D) in enumerate([
                ("Aitchison", _D_aitchison),
                ("Jensen-Shannon", _D_js),
            ]):
                _sub = _D[np.ix_(
                    _sample_indices,
                    _sample_indices,
                )]
                _D2 = _sub ** 2

                # Different but reproducible seeds for each test
                _seed = 10_000 * _distance_index + _var_index

                _F, _R2, _p_perm = permanova(
                    _D2,
                    _codes.copy(),
                    perms=perms,
                    seed=_seed,
                )

                _dispersion_F, _dispersion_p, _ = permdisp(
                    _D2,
                    _codes.copy(),
                    perms=perms,
                    seed=_seed,
                )

                _rows.append({
                    "variable": _var,
                    "distance": _name,
                    "n": int(_keep_mask.sum()),
                    "groups": int(len(np.unique(_codes))),
                    "pseudo_F": _F,
                    "R2": _R2,
                    "p_perm": _p_perm,
                    "dispersion_F": _dispersion_F,
                    "dispersion_p": _dispersion_p,
                })

        return (
            pd.DataFrame(_rows)
            .sort_values(
                ["distance", "R2"],
                ascending=[True, False],
            )
            .reset_index(drop=True)
        )

    permanova_results = _run_permanova()

    permanova_results.round(4)
    return (permanova_results,)


@app.cell
def _(meta_vars, mo):

    # Interactive variable selector

    var_pick = mo.ui.dropdown(
        options=meta_vars,
        value=meta_vars[0] if meta_vars else None,
        label="Variable:",
    )

    var_pick
    return (var_pick,)


@app.cell
def _(P, meta, np, plt, re, topic_cols, topic_labels, var_pick):

    # Mean microbial-component composition by group

    def _composition():
        _variable = var_pick.value

        if _variable is None:
            _fig, _ax = plt.subplots(figsize=(6, 2))
            _ax.text(0.5, 0.5, "No testable variable", ha="center")
            _ax.axis("off")
            return _fig

        _labels = meta[_variable]

        _categories = [
            category
            for category in sorted(
                _labels[_labels != "n/a"].unique()
            )
            if (_labels == category).sum() >= 3
        ]

        _sample_sizes = [
            int((_labels == category).sum())
            for category in _categories
        ]

        def _shorten(items):
            _base = [
                re.sub(r"\s*\(.*?\)", "", item).strip()
                for item in items
            ]

            if len(_base) > 1:
                _tokens = [item.split(":") for item in _base]
                _minimum_length = min(len(tokens) for tokens in _tokens)
                _common = 0

                for _i in range(max(0, _minimum_length - 1)):
                    _segment = _tokens[0][_i].strip()
                    if all(
                        tokens[_i].strip() == _segment
                        for tokens in _tokens
                    ):
                        _common += 1
                    else:
                        break

                if _common:
                    _base = [
                        ":".join(tokens[_common:])
                        for tokens in _tokens
                    ]

            _output = []

            for item in _base:
                _clean = re.sub(r"\s+", " ", item).strip(" :")
                _output.append(
                    (_clean[:25] + "…")
                    if len(_clean) > 26
                    else (_clean or "?")
                )

            return _output

        _display_labels = [
            f"{short_label}  (n={sample_size})"
            for short_label, sample_size in zip(
                _shorten(_categories),
                _sample_sizes,
            )
        ]

        # Mean topic proportion within each group
        _mean_composition = np.zeros(
            (len(_categories), len(topic_cols)),
            dtype=float,
        )

        for _i, _category in enumerate(_categories):
            _mean = (
                P.loc[
                    (_labels == _category).to_numpy(),
                    topic_cols,
                ]
                .mean(axis=0)
                .to_numpy()
            )

            _total = _mean.sum()
            _mean_composition[_i] = (
                _mean / _total if _total > 0 else _mean
            )

        _colours = [
            "#4C72B0",
            "#DD8452",
            "#55A868",
            "#C44E52",
            "#8172B3",
            "#937860",
            "#DA8BC3",
            "#8C8C8C",
            "#CCB974",
            "#64B5CD",
        ]

        _fig, _ax = plt.subplots(
            figsize=(11, 0.62 * len(_categories) + 1.8),
            constrained_layout=True,
        )

        _y_positions = np.arange(len(_categories))
        _left = np.zeros(len(_categories))

        for _topic_index in range(len(topic_cols)):
            _ax.barh(
                _y_positions,
                _mean_composition[:, _topic_index],
                left=_left,
                height=0.72,
                color=_colours[_topic_index % len(_colours)],
                edgecolor="white",
                linewidth=0.6,
                label=topic_labels[_topic_index],
                zorder=3,
            )

            _left += _mean_composition[:, _topic_index]

        _ax.set_yticks(_y_positions)
        _ax.set_yticklabels(_display_labels, fontsize=8)
        _ax.invert_yaxis()
        _ax.set_xlim(0, 1)
        _ax.set_xlabel("Mean microbial-component proportion")
        _ax.set_title(
            f"Mean microbial-component composition by {_variable}"
        )

        _ax.legend(
            loc="center left",
            bbox_to_anchor=(1.01, 0.5),
            frameon=False,
            fontsize=8,
            title="Component",
        )

        return _fig

    _composition()
    return


@app.cell
def _(P, clr, jensenshannon, np):

    # Distance matrices shared by the interactive plots

    def _calculate_distances():
        from scipy.spatial.distance import pdist
        from scipy.spatial.distance import squareform

        _P_normalised = P.div(P.sum(axis=1), axis=0).to_numpy()

        _D_aitchison = squareform(
            pdist(
                clr(_P_normalised),
                metric="euclidean",
            )
        )

        _n = len(_P_normalised)
        _D_js = np.zeros((_n, _n), dtype=float)

        for _i in range(_n):
            for _j in range(_i + 1, _n):
                _distance = jensenshannon(
                    _P_normalised[_i],
                    _P_normalised[_j],
                )

                _distance = (
                    0.0 if np.isnan(_distance) else float(_distance)
                )

                _D_js[_i, _j] = _distance
                _D_js[_j, _i] = _distance

        return {
            "Aitchison": _D_aitchison,
            "Jensen-Shannon": _D_js,
        }

    DIST = _calculate_distances()
    return (DIST,)


@app.cell
def _(mo):

    # Interactive distance selector

    dist_pick = mo.ui.dropdown(
        options=["Aitchison", "Jensen-Shannon"],
        value="Jensen-Shannon",
        label="Distance:",
    )

    dist_pick
    return (dist_pick,)


@app.cell
def _(DIST, dist_pick, meta, np, pd, permanova, permdisp, plt, re, var_pick):

    # PCoA ordination and dispersion plot

    def _ordination():
        from matplotlib.patches import Ellipse

        _variable = var_pick.value
        _distance_name = dist_pick.value

        if _variable is None:
            _fig, _ax = plt.subplots(figsize=(6, 2))
            _ax.text(0.5, 0.5, "No testable variable", ha="center")
            _ax.axis("off")
            return _fig

        _D = DIST[_distance_name]
        _labels = meta[_variable].to_numpy()

        _nonmissing = _labels != "n/a"
        _counts = pd.Series(
            _labels[_nonmissing]
        ).value_counts()

        _keep = _counts[_counts >= 3].index
        _keep_mask = (
            _nonmissing
            & pd.Series(_labels).isin(_keep).to_numpy()
        )

        _sample_indices = np.where(_keep_mask)[0]
        _category = pd.Categorical(_labels[_keep_mask])
        _categories = list(_category.categories)

        if len(_categories) < 2:
            _fig, _ax = plt.subplots(figsize=(6, 2))
            _ax.text(0.5, 0.5, "Need at least two groups", ha="center")
            _ax.axis("off")
            return _fig

        _sub = _D[np.ix_(
            _sample_indices,
            _sample_indices,
        )]
        _D2 = _sub ** 2

        _n = len(_sample_indices)
        _J = np.eye(_n) - np.ones((_n, _n)) / _n
        _B = -0.5 * _J @ _D2 @ _J

        _eigenvalues, _eigenvectors = np.linalg.eigh(_B)
        _order = np.argsort(_eigenvalues)[::-1]
        _eigenvalues = _eigenvalues[_order]
        _eigenvectors = _eigenvectors[:, _order]

        _positive = _eigenvalues > 1e-10
        _positive_values = _eigenvalues[_positive]
        _coordinates = (
            _eigenvectors[:, _positive]
            * np.sqrt(_positive_values)
        )

        if _coordinates.shape[1] == 1:
            _coordinates = np.column_stack([
                _coordinates[:, 0],
                np.zeros(_coordinates.shape[0]),
            ])
            _axis_fraction = np.array([1.0, 0.0])
        else:
            _axis_fraction = (
                _positive_values / _positive_values.sum()
            )

        _XY = _coordinates[:, :2]

        _F, _R2, _p_perm = permanova(
            _D2,
            _category.codes.copy(),
            perms=9999,
            seed=0,
        )

        _dispersion_F, _dispersion_p, _distance_to_centroid = permdisp(
            _D2,
            _category.codes.copy(),
            perms=9999,
            seed=0,
        )

        def _shorten(items):
            _base = [
                re.sub(r"\s*\(.*?\)", "", item).strip()
                for item in items
            ]

            if len(_base) > 1:
                _tokens = [item.split(":") for item in _base]
                _minimum_length = min(
                    len(tokens) for tokens in _tokens
                )
                _common = 0

                for _i in range(max(0, _minimum_length - 1)):
                    if all(
                        tokens[_i].strip()
                        == _tokens[0][_i].strip()
                        for tokens in _tokens
                    ):
                        _common += 1
                    else:
                        break

                if _common:
                    _base = [
                        ":".join(tokens[_common:])
                        for tokens in _tokens
                    ]

            return [
                re.sub(r"\s+", " ", item).strip(" :")[:24]
                for item in _base
            ]

        _short_labels = _shorten(_categories)

        _palette = [
            "#4C72B0",
            "#DD8452",
            "#55A868",
            "#C44E52",
            "#8172B3",
            "#937860",
            "#DA8BC3",
            "#8C8C8C",
            "#CCB974",
            "#64B5CD",
        ]

        _fig, (_axis_ordination, _axis_dispersion) = plt.subplots(
            1,
            2,
            figsize=(12.5, 5),
            gridspec_kw={"width_ratios": [1.5, 1]},
        )

        for _group_index in range(len(_categories)):
            _group_rows = np.where(
                _category.codes == _group_index
            )[0]

            _colour = _palette[
                _group_index % len(_palette)
            ]

            _axis_ordination.scatter(
                _XY[_group_rows, 0],
                _XY[_group_rows, 1],
                s=18,
                color=_colour,
                alpha=0.6,
                edgecolors="none",
                label=(
                    f"{_short_labels[_group_index]} "
                    f"(n={len(_group_rows)})"
                ),
            )

            _centroid = _XY[_group_rows].mean(axis=0)

            _axis_ordination.scatter(
                *_centroid,
                color=_colour,
                s=130,
                marker="X",
                edgecolors="black",
                linewidths=0.8,
                zorder=5,
            )

            if len(_group_rows) > 2:
                _covariance = np.cov(
                    _XY[_group_rows].T
                )

                _ellipse_values, _ellipse_vectors = np.linalg.eigh(
                    _covariance
                )
                _ellipse_order = _ellipse_values.argsort()[::-1]
                _ellipse_values = _ellipse_values[_ellipse_order]
                _ellipse_vectors = _ellipse_vectors[:, _ellipse_order]

                _angle = np.degrees(
                    np.arctan2(
                        _ellipse_vectors[1, 0],
                        _ellipse_vectors[0, 0],
                    )
                )

                _ellipse = Ellipse(
                    _centroid,
                    4 * np.sqrt(max(_ellipse_values[0], 0)),
                    4 * np.sqrt(max(_ellipse_values[1], 0)),
                    angle=_angle,
                    facecolor=_colour,
                    alpha=0.12,
                    edgecolor=_colour,
                    lw=1.2,
                )

                _axis_ordination.add_patch(_ellipse)

        _axis_ordination.axhline(
            0,
            color="#dddddd",
            lw=0.6,
        )
        _axis_ordination.axvline(
            0,
            color="#dddddd",
            lw=0.6,
        )

        _axis_ordination.set_xlabel(
            f"PCoA1 ({_axis_fraction[0] * 100:.1f}%)"
        )
        _axis_ordination.set_ylabel(
            f"PCoA2 ({_axis_fraction[1] * 100:.1f}%)"
        )

        _axis_ordination.set_title(
            f"{_variable} · {_distance_name}\n"
            f"PERMANOVA R²={_R2:.3f}, p={_p_perm:.4f}"
        )

        _axis_ordination.legend(
            fontsize=7,
            loc="best",
            frameon=False,
        )

        _rng = np.random.default_rng(0)

        for _group_index in range(len(_categories)):
            _group_rows = np.where(
                _category.codes == _group_index
            )[0]

            _colour = _palette[
                _group_index % len(_palette)
            ]

            _jitter = (
                _rng.random(len(_group_rows)) - 0.5
            ) * 0.3

            _axis_dispersion.scatter(
                _distance_to_centroid[_group_rows],
                np.full(
                    len(_group_rows),
                    _group_index,
                ) + _jitter,
                s=12,
                color=_colour,
                alpha=0.5,
                edgecolors="none",
            )

            _median = np.median(
                _distance_to_centroid[_group_rows]
            )

            _axis_dispersion.plot(
                [_median, _median],
                [
                    _group_index - 0.28,
                    _group_index + 0.28,
                ],
                color="black",
                lw=1.6,
            )

        _axis_dispersion.set_yticks(
            range(len(_categories))
        )
        _axis_dispersion.set_yticklabels(
            _short_labels,
            fontsize=7,
        )
        _axis_dispersion.set_ylim(
            len(_categories) - 0.5,
            -0.5,
        )
        _axis_dispersion.set_xlabel(
            "Distance to group centroid"
        )
        _axis_dispersion.set_title(
            "PERMDISP\n"
            f"F={_dispersion_F:.3f}, p={_dispersion_p:.4f}"
        )

        _fig.tight_layout()
        return _fig

    _ordination()
    return


@app.cell
def _(dist_pick, np, permanova_results, plt, significance_stars):

    # PERMANOVA R² summary

    def _r2_bar():
        _distance_name = dist_pick.value

        _table = (
            permanova_results[
                permanova_results["distance"] == _distance_name
            ]
            .sort_values("R2")
            .reset_index(drop=True)
        )

        _fig, _ax = plt.subplots(
            figsize=(7.7, 0.52 * len(_table) + 1.5)
        )

        _y = np.arange(len(_table))

        for _i, _row in _table.iterrows():
            _ax.barh(
            _i,
            _row["R2"],
            color="#4C72B0",
            alpha=0.85,
            edgecolor="white",
            zorder=3,
        )

        _star = significance_stars(_row["p_perm"])
        _label = (
            f'{_row["R2"]:.3f} {_star}'
            if _star
            else f'{_row["R2"]:.3f}'
        )
        _ax.text(_row["R2"] + 0.003, _i, _label, va="center", fontsize=8)

        _ax.set_yticks(_y)
        _ax.set_yticklabels(
            _table["variable"],
            fontsize=9,
        )

        _ax.set_xlabel(
            "R² (share of compositional variance)"
        )

        _ax.set_title(
            f"PERMANOVA R² by variable · {_distance_name}\n"
            "* p≤.05  ** p≤.01  *** p≤.001"
            "  ·  hatched = significant PERMDISP"
            " among significant PERMANOVA effects"
        )

        _maximum = (
            _table["R2"].max()
            if not _table.empty
            else 0.05
        )

        _ax.set_xlim(
            0,
            max(0.05, _maximum * 1.30),
        )

        _fig.tight_layout()
        return _fig

    _r2_bar()
    return


@app.cell
def _():
    return


@app.cell
def _(DIST, dist_pick, meta, np, pd, permanova, permdisp, plt, re, var_pick):
    # PCoA ordination and dispersion plot
    # Only the plotting style has changed; the PERMANOVA / PERMDISP calls,
    # distances, group coding and permutation seeds are identical.

    def _ordination():
        from matplotlib.patches import Ellipse

        _variable = var_pick.value
        _distance_name = dist_pick.value

        if _variable is None:
            _fig, _ax = plt.subplots(figsize=(6, 2))
            _ax.text(0.5, 0.5, "No testable variable", ha="center")
            _ax.axis("off")
            return _fig

        _D = DIST[_distance_name]
        _labels = meta[_variable].to_numpy()

        _nonmissing = _labels != "n/a"
        _counts = pd.Series(
            _labels[_nonmissing]
        ).value_counts()

        _keep = _counts[_counts >= 3].index
        _keep_mask = (
            _nonmissing
            & pd.Series(_labels).isin(_keep).to_numpy()
        )

        _sample_indices = np.where(_keep_mask)[0]
        _category = pd.Categorical(_labels[_keep_mask])
        _categories = list(_category.categories)

        if len(_categories) < 2:
            _fig, _ax = plt.subplots(figsize=(6, 2))
            _ax.text(0.5, 0.5, "Need at least two groups", ha="center")
            _ax.axis("off")
            return _fig

        _sub = _D[np.ix_(
            _sample_indices,
            _sample_indices,
        )]
        _D2 = _sub ** 2

        _n = len(_sample_indices)
        _J = np.eye(_n) - np.ones((_n, _n)) / _n
        _B = -0.5 * _J @ _D2 @ _J

        _eigenvalues, _eigenvectors = np.linalg.eigh(_B)
        _order = np.argsort(_eigenvalues)[::-1]
        _eigenvalues = _eigenvalues[_order]
        _eigenvectors = _eigenvectors[:, _order]

        _positive = _eigenvalues > 1e-10
        _positive_values = _eigenvalues[_positive]
        _coordinates = (
            _eigenvectors[:, _positive]
            * np.sqrt(_positive_values)
        )

        if _coordinates.shape[1] == 1:
            _coordinates = np.column_stack([
                _coordinates[:, 0],
                np.zeros(_coordinates.shape[0]),
            ])
            _axis_fraction = np.array([1.0, 0.0])
        else:
            _axis_fraction = (
                _positive_values / _positive_values.sum()
            )

        _XY = _coordinates[:, :2]

        _F, _R2, _p_perm = permanova(
            _D2,
            _category.codes.copy(),
            perms=9999,
            seed=0,
        )

        _dispersion_F, _dispersion_p, _distance_to_centroid = permdisp(
            _D2,
            _category.codes.copy(),
            perms=9999,
            seed=0,
        )

        def _shorten(items):
            _base = [
                re.sub(r"\s*\(.*?\)", "", item).strip()
                for item in items
            ]

            if len(_base) > 1:
                _tokens = [item.split(":") for item in _base]
                _minimum_length = min(
                    len(tokens) for tokens in _tokens
                )
                _common = 0

                for _i in range(max(0, _minimum_length - 1)):
                    if all(
                        tokens[_i].strip()
                        == _tokens[0][_i].strip()
                        for tokens in _tokens
                    ):
                        _common += 1
                    else:
                        break

                if _common:
                    _base = [
                        ":".join(tokens[_common:])
                        for tokens in _tokens
                    ]

            return [
                re.sub(r"\s+", " ", item).strip(" :")[:24]
                for item in _base
            ]

        _short_labels = _shorten(_categories)

        # Display names only (legend and PERMDISP axis). The underlying
        # group labels used for the statistics are not changed.
        def _display_name(label):
            _l = label.lower()
            if "pre" in _l and "expos" in _l:
                return "Pre-exposure"
            if "long" in _l or "survivor" in _l:
                return "Long-term survivors"
            if "juvenile" in _l:
                return "23-day exposure – juvenile"
            if "expos" in _l or "post" in _l:
                return "23-day exposure – adult"
            return label

        _short_labels = [_display_name(label) for label in _short_labels]

        # Display name for the variable in the plot title
        _variable_title = {"Category": "Exposure time"}.get(_variable, _variable)

        # Colour scheme: pre-exposure = green, 23-day exposure adult = red,
        # juveniles = blue, long-term survivors = orange.
        _fallback_palette = [
            "#8172B3",
            "#937860",
            "#DA8BC3",
            "#8C8C8C",
            "#CCB974",
            "#64B5CD",
        ]

        def _group_colour(label, index):
            _l = label.lower()
            if "pre" in _l and "expos" in _l:
                return "#55A868"   # green
            if "survivor" in _l or "long" in _l:
                return "#E0A458"   # muted orange (more yellow, so it stays distinct from red)
            if "juvenile" in _l:
                return "#4C72B0"   # blue
            if "expos" in _l or "post" in _l:
                return "#C44E52"   # red
            return _fallback_palette[index % len(_fallback_palette)]

        _colours = [
            _group_colour(label, i)
            for i, label in enumerate(_categories)
        ]

        _fig, (_axis_ordination, _axis_dispersion) = plt.subplots(
            1,
            2,
            figsize=(12.5, 5),
            gridspec_kw={"width_ratios": [1.5, 1]},
        )

        for _group_index in range(len(_categories)):
            _group_rows = np.where(
                _category.codes == _group_index
            )[0]

            _colour = _colours[_group_index]

            _axis_ordination.scatter(
                _XY[_group_rows, 0],
                _XY[_group_rows, 1],
                s=18,
                color=_colour,
                alpha=0.6,
                edgecolors="none",
                label=(
                    f"{_short_labels[_group_index]} "
                    f"(n={len(_group_rows)})"
                ),
            )

            _centroid = _XY[_group_rows].mean(axis=0)

            _axis_ordination.scatter(
                *_centroid,
                color=_colour,
                s=130,
                marker="X",
                edgecolors="black",
                linewidths=0.8,
                zorder=5,
            )

            if len(_group_rows) > 2:
                _covariance = np.cov(
                    _XY[_group_rows].T
                )

                _ellipse_values, _ellipse_vectors = np.linalg.eigh(
                    _covariance
                )
                _ellipse_order = _ellipse_values.argsort()[::-1]
                _ellipse_values = _ellipse_values[_ellipse_order]
                _ellipse_vectors = _ellipse_vectors[:, _ellipse_order]

                _angle = np.degrees(
                    np.arctan2(
                        _ellipse_vectors[1, 0],
                        _ellipse_vectors[0, 0],
                    )
                )

                _ellipse = Ellipse(
                    _centroid,
                    4 * np.sqrt(max(_ellipse_values[0], 0)),
                    4 * np.sqrt(max(_ellipse_values[1], 0)),
                    angle=_angle,
                    facecolor=_colour,
                    alpha=0.12,
                    edgecolor=_colour,
                    lw=1.2,
                )

                _axis_ordination.add_patch(_ellipse)

        _axis_ordination.axhline(
            0,
            color="#dddddd",
            lw=0.6,
        )
        _axis_ordination.axvline(
            0,
            color="#dddddd",
            lw=0.6,
        )

        _axis_ordination.set_xlabel(
            f"PCoA1 ({_axis_fraction[0] * 100:.1f}%)"
        )
        _axis_ordination.set_ylabel(
            f"PCoA2 ({_axis_fraction[1] * 100:.1f}%)"
        )

        _axis_ordination.set_title(
            f"{_variable_title} · {_distance_name}\n"
            f"PERMANOVA R²={_R2:.3f}, p={_p_perm:.4f}"
        )

        _axis_ordination.legend(
            fontsize=7,
            loc="best",
            frameon=False,
        )

        # ---------------------------------------------------------------
        # PERMDISP panel, laid out in two blocks:
        #   Juvenile | 23-day exposure
        #   ---------------------------------------
        #   Adult    | Long-term survivors
        #            | 23-day exposure
        #            | Pre-exposure
        # Only the row order and labels change; distances are unchanged.
        # ---------------------------------------------------------------
        def _row_info(raw_label):
            _l = raw_label.lower()
            if "juvenile" in _l:
                return 0, "Juvenile", "23-day exposure"
            if "long" in _l or "survivor" in _l:
                return 1, "Adult", "Long-term survivors"
            if "pre" in _l and "expos" in _l:
                return 3, "Adult", "Pre-exposure"
            if "expos" in _l or "post" in _l:
                return 2, "Adult", "23-day exposure"
            return 4, "Other", raw_label

        _info = [_row_info(label) for label in _categories]
        _sorted_groups = sorted(
            range(len(_categories)), key=lambda g: (_info[g][0], g)
        )
        _row_of = {g: row for row, g in enumerate(_sorted_groups)}

        _rng = np.random.default_rng(0)

        for _group_index in range(len(_categories)):
            _group_rows = np.where(
                _category.codes == _group_index
            )[0]

            _colour = _colours[_group_index]
            _y = _row_of[_group_index]

            _jitter = (
                _rng.random(len(_group_rows)) - 0.5
            ) * 0.3

            _axis_dispersion.scatter(
                _distance_to_centroid[_group_rows],
                np.full(len(_group_rows), _y) + _jitter,
                s=12,
                color=_colour,
                alpha=0.5,
                edgecolors="none",
            )

            _median = np.median(
                _distance_to_centroid[_group_rows]
            )

            _axis_dispersion.plot(
                [_median, _median],
                [_y - 0.28, _y + 0.28],
                color="black",
                lw=1.6,
            )

        _n_rows = len(_categories)
        _axis_dispersion.set_yticks(range(_n_rows))
        _axis_dispersion.set_yticklabels(
            [_info[g][2] for g in _sorted_groups],
            fontsize=8,
        )
        _axis_dispersion.set_ylim(_n_rows - 0.5, -0.5)
        _axis_dispersion.set_xlabel(
            "Distance to group centroid"
        )
        _axis_dispersion.set_title(
            "PERMDISP\n"
            f"F={_dispersion_F:.3f}, p={_dispersion_p:.4f}"
        )

        # Stage labels ("Juvenile", "Adult") sit left of the tick labels.
        # Offsets are in points, so they stay put when the layout changes.
        from matplotlib.lines import Line2D
        from matplotlib.transforms import (
            blended_transform_factory,
            offset_copy,
        )

        _fig.canvas.draw()
        _renderer = _fig.canvas.get_renderer()
        _px_to_pt = 72.0 / _fig.dpi
        _tick_width_pt = _px_to_pt * max(
            t.get_window_extent(_renderer).width
            for t in _axis_dispersion.get_yticklabels()
        )
        _tick_pad_pt = 3.5 + 2.0  # default tick length + pad

        _base = blended_transform_factory(
            _axis_dispersion.transAxes, _axis_dispersion.transData
        )
        _bracket_pt = _tick_pad_pt + _tick_width_pt + 6
        _stage_pt = _bracket_pt + 5
        _trans_bracket = offset_copy(
            _base, fig=_fig, x=-_bracket_pt, y=0, units="points"
        )
        _trans_stage = offset_copy(
            _base, fig=_fig, x=-_stage_pt, y=0, units="points"
        )

        # Stage blocks: consecutive rows sharing a stage
        _stages = [_info[g][1] for g in _sorted_groups]
        _blocks = []
        for _row, _stage in enumerate(_stages):
            if _blocks and _blocks[-1][0] == _stage:
                _blocks[-1][2] = _row
            else:
                _blocks.append([_stage, _row, _row])

        _stage_texts = []
        for _stage, _first, _last in _blocks:
            _stage_texts.append(_axis_dispersion.text(
                0,
                (_first + _last) / 2,
                _stage,
                transform=_trans_stage,
                ha="right",
                va="center",
                fontsize=9,
                clip_on=False,
            ))
            if _last > _first:
                _axis_dispersion.plot(
                    [0, 0],
                    [_first - 0.35, _last + 0.35],
                    transform=_trans_bracket,
                    color="black",
                    lw=0.8,
                    clip_on=False,
                )

        # Separator line between blocks, running under the labels too
        _fig.canvas.draw()
        _stage_width_pt = _px_to_pt * max(
            t.get_window_extent(_renderer).width for t in _stage_texts
        )
        _left_pt = _stage_pt + _stage_width_pt + 4
        for _stage, _first, _last in _blocks[1:]:
            _y_sep = _first - 0.5
            # across the plot ...
            _axis_dispersion.add_line(Line2D(
                [0, 1], [_y_sep, _y_sep],
                transform=_base, color="black", lw=0.8, clip_on=False,
            ))
            # ... and continued left under the labels
            _axis_dispersion.annotate(
                "",
                xy=(0, _y_sep), xycoords=_base,
                xytext=(-_left_pt, 0), textcoords="offset points",
                arrowprops=dict(arrowstyle="-", color="black", lw=0.8,
                                shrinkA=0, shrinkB=0),
                annotation_clip=False,
            )

        _fig.tight_layout()

        # ---------------------------------------------------------------
        # Export: PDF and SVG (vector, editable text) plus a 300 dpi PNG.
        # Files are written next to the notebook.
        # ---------------------------------------------------------------
        _stem = re.sub(
            r"[^A-Za-z0-9]+",
            "_",
            f"PCoA_PERMDISP_{_variable}_{_distance_name}",
        ).strip("_")

        with plt.rc_context({"pdf.fonttype": 42, "svg.fonttype": "none"}):
            for _ext in ("pdf", "svg", "png"):
                _fig.savefig(
                    f"{_stem}.{_ext}",
                    dpi=300,
                    bbox_inches="tight",
                )

        return _fig

    _ordination()
    return


if __name__ == "__main__":
    app.run()
