def L_scv_pl_heatmap(
    adata,
    var_names,
    sortby="latent_time",
    layer="Ms",
    color_map="viridis",
    col_color=None,
    palette="viridis",
    n_convolve=30,
    standard_scale=0,
    sort=True,
    colorbar=None,
    col_cluster=False,
    row_cluster=False,
    context=None,
    font_scale=None,
    figsize=(8, 4),
    show=None,
    save=None,
    debug=False,
    label_genes=None,
    min_label_spacing = 0.1,
    x_outside = 1.1,
    print_actual_gene_order = False,
    **kwargs,
):
    """
    Plot a heatmap of gene expression across latent time using seaborn's clustermap,
    with optional smoothing, sorting, column color annotations, and external gene labels.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix from scVelo or Scanpy with expression and metadata.
    
    var_names : list of str
        List of gene names to include in the heatmap (must exist in `adata.var_names`).
    
    sortby : str, default "latent_time"
        Observation key to use for sorting cells (e.g., pseudotime).
    
    layer : str, default "Ms"
        Layer of `adata` to use for expression values (e.g., spliced, unspliced, total).
    
    color_map : str, default "viridis"
        Colormap to use for the heatmap.
    
    col_color : str or list of str, optional
        Keys from `adata.obs` to add as categorical top bars in the heatmap.
    
    palette : str, default "viridis"
        Color palette to use for categorical color annotations.
    
    n_convolve : int or None, default 30
        If set, smooth expression along pseudotime using a moving average of this window size.
    
    standard_scale : int (0 or 1), default 0
        Whether to standardize across rows (0) or columns (1).
    
    sort : bool, default True
        Whether to sort genes by the peak time of expression.
    
    colorbar : bool or None
        Whether to display the expression colorbar.
    
    col_cluster : bool, default False
        Whether to cluster columns (cells).
    
    row_cluster : bool, default False
        Whether to cluster rows (genes).
    
    context : str or None
        Seaborn context (e.g., "notebook", "talk", "poster") for font scaling.
    
    font_scale : float or None
        Font scaling factor for seaborn context.
    
    figsize : tuple, default (8, 4)
        Size of the figure in inches.
    
    show : bool or None
        Whether to show the plot immediately.
    
    save : str or None
        Path to save the plot if desired.
    
    debug : bool, default False
        If True, prints intermediate data for debugging.
    
    label_genes : list of str or None
        Subset of genes to label outside the heatmap with connecting lines.
    
    min_label_spacing : float, default 0.1
        Minimum vertical spacing between gene labels (in Axes coordinates).
    
    x_outside : float, default 1.1
        Horizontal position to place gene labels outside the heatmap (Axes coordinates).
    
    print_actual_gene_order : bool, default False
        If True, prints the actual gene order used in the heatmap after clustering.

    Returns
    -------
    cm : seaborn.matrix.ClusterGrid or None
        The ClusterGrid object representing the heatmap (only returned if show=False).
    """
        
    import seaborn as sns
    import numpy as np
    import pandas as pd
    from scipy.sparse import issparse

    from scvelo import logging as logg
    from scvelo.tools.utils import strings_to_categoricals
    from scvelo.utils import (
        interpret_colorkey,
        is_categorical,
    )
    from scvelo.plotting.utils import (
        savefig_or_show,
        set_colors_for_categorical_obs,
        to_list,
    )

    var_names = [name for name in var_names if name in adata.var_names]
    if not var_names:
        logg.error("No valid genes found in `adata.var_names` for plotting. Please check your `var_names` input.")
        return # Or raise an error

    tkey, xkey = kwargs.pop("tkey", sortby), kwargs.pop("xkey", layer)
    time = adata.obs[tkey].values
    valid = np.isfinite(time)
    if not np.any(valid):
        logg.error(f"No finite values found in adata.obs['{tkey}']. Cannot sort and plot.")
        return
    time = time[valid]

    X = (
        adata[:, var_names].layers[xkey]
        if xkey in adata.layers.keys()
        else adata[:, var_names].X
    )
    if issparse(X):
        X = X.toarray()
    X = X[valid]
    if debug:
        print("X is:")
        print(X)
    df = pd.DataFrame(X[np.argsort(time)], columns=var_names)
    if debug:
        print('before setting index')
        print(df)
    df.index = adata.obs_names[valid][np.argsort(time)]
    if debug:
        print('after setting index')
        print(df)

    if n_convolve is not None:
        weights = np.ones(n_convolve) / n_convolve
        for gene in var_names:
            try:
                df[gene] = np.convolve(df[gene].values, weights, mode="same")
            except ValueError as e:
                logg.info(f"Skipping variable {gene}: {e}")
                pass

    if sort:
        max_sort = np.argsort(np.argmax(df.values, axis=0))
        if debug:
            print('if sort w/o index')
            print(pd.DataFrame(df.values[:, max_sort], columns=df.columns[max_sort]))
        df = pd.DataFrame(df.values[:, max_sort], columns=df.columns[max_sort], index=df.index)
        if debug:
            print('sorted w index')
            print(df)

    strings_to_categoricals(adata)

    #<===== Handle multiple col_color keys as a DataFrame of color tracks =====>
    if col_color is not None:
        original_col_color_keys = to_list(col_color)
        col_color_dict = {}
    else:
        original_col_color_keys = None

    if original_col_color_keys is not None:
        for col in original_col_color_keys:
            if not is_categorical(adata, col):
                obs_col = adata.obs[col]
                cat_col = np.round(obs_col / np.max(obs_col), 2) * np.max(obs_col)
                adata.obs[f"{col}_categorical"] = pd.Categorical(cat_col)
                col += "_categorical"
                set_colors_for_categorical_obs(adata, col, palette)
            col_color_dict[col] = interpret_colorkey(adata, col)[np.argsort(time)]
        col_colors = pd.DataFrame(col_color_dict, index=adata.obs_names[valid][np.argsort(time)])# Combine all into a DataFrame for seaborn
    
    #<=== below set cbar_pos required for colorbar ===>
    if "cbar_pos" not in kwargs:
        kwargs["cbar_pos"] = (0.02, 0.8, 0.03, 0.18) if colorbar is not False else None
    kwargs.setdefault("cbar_kws", {"label": "expression"})

    args = {}
    if font_scale is not None:
        args = {"font_scale": font_scale}
        context = context or "notebook"

    with sns.plotting_context(context=context, **args):
        if debug:
            print(f"Shape of col_colors DataFrame: {col_colors.shape}")
            print(f"Number of NaN values in col_colors: {col_colors.isnull().sum().sum()}")
            print(f"Are col_colors empty? {col_colors.empty}")
            print(f"col_colors: {col_colors}")
            print(f"DataFrame shape before clustermap: {df.shape}")
            print(f"Number of NaN values in df: {df.isnull().sum().sum()}")
        
        #<=== below enable col_colors is rendered as the top bar ===>
        col_colors = col_colors.loc[df.index]
        
        #<=== remove the label of catagory bar "col_colors" ===>
        col_colors.columns = [""]
        
        kwargs.update(
            {
                "col_colors": col_colors,
                "col_cluster": col_cluster,
                "row_cluster": row_cluster,
                "cmap": color_map,
                "xticklabels": False,
                "standard_scale": standard_scale,
                "figsize": figsize,
            }
        )

        try:
            cm = sns.clustermap(df.T, **kwargs)
        except ImportWarning:
            logg.warn("Please upgrade seaborn with `pip install -U seaborn`.")
            kwargs.pop("dendrogram_ratio")
            kwargs.pop("cbar_pos")
            cm = sns.clustermap(df.T, **kwargs)
        
        #<=== remove top tick marks of catagory bar "col_colors" ===>
        if hasattr(cm, "ax_col_colors"):
            cm.ax_col_colors.set_yticks([])
        
        #<=== Specify the genes to label ===>
        L_add_gene_labels(cm, df, label_genes,
                           min_label_spacing = min_label_spacing,
                           x_outside = x_outside,
                           print_actual_gene_order = print_actual_gene_order
                          )
       
        #<=== Manually add legends for categorical col_colors ===>
        if original_col_color_keys is not None:
            import matplotlib.patches as mpatches
            for col, color_array in zip(original_col_color_keys, col_colors):
                if is_categorical(adata, col):
                    categories = adata.obs[col].cat.categories
                    colors = adata.uns.get(f"{col}_colors", [])
                    if categories is not None and colors is not None:
                        handles = [
                            mpatches.Patch(color=c, label=cat)
                            for cat, c in zip(categories, colors)
                        ]
                        cm.ax_col_dendrogram.legend(
                            handles=handles,
                            title=col,
                            bbox_to_anchor=(1.05, 1),
                            loc="upper left",
                            fontsize="small",
                            frameon=False,
                        )

    savefig_or_show("heatmap", save=save, show=show, dpi=300)
    if show is False:
        return cm


def L_add_gene_labels(cm, df, label_genes,
                        min_label_spacing,
                        x_outside,
                        print_actual_gene_order
                          ):
    """
    Add gene labels outside the heatmap along the Y-axis with repelling to reduce overlap.

    This function uses Axes coordinates to place gene names outside the heatmap and
    draws connecting lines to the actual rows. It avoids overlap by pushing labels apart
    vertically while preserving order.

    Parameters
    ----------
    cm : seaborn.matrix.ClusterGrid
        The heatmap object generated by seaborn.clustermap.

    df : pd.DataFrame
        The gene expression data used in the heatmap, with genes in columns.

    label_genes : list of str
        Genes to label outside the heatmap. Only those in `df.columns` are considered.

    min_label_spacing : float
        Minimum spacing between labels (in Axes coordinate units). Helps avoid overlap.

    x_outside : float
        Horizontal position to place labels (Axes coordinate system, >1 means right of heatmap).

    print_actual_gene_order : bool
        If True, prints the actual row order of genes after any reordering done by clustering.

    Returns
    -------
    None
        This function modifies the heatmap in-place by adding text and annotation lines.
    """
    if label_genes is None:
        return
        
    label_set = set(label_genes)
    
    # Get the actual row order from the clustermap object
    # This accounts for any reordering done by seaborn during plotting
    if hasattr(cm, 'dendrogram_row') and cm.dendrogram_row is not None:
        # If row clustering was used, get the reordered indices
        row_order = cm.dendrogram_row.reordered_ind
        actual_gene_order = [df.columns[i] for i in row_order]
    else:
        # No row clustering, use the original order
        actual_gene_order = df.columns.to_list()
    
    if print_actual_gene_order:
        print(actual_gene_order)
    # Find positions of labeled genes in the ACTUAL display order
    labeled_genes_data = []
    for i, gene in enumerate(actual_gene_order):
        if gene in label_set:
            gene_y_coord = 1.0 - (i + 0.5) / len(actual_gene_order)
            labeled_genes_data.append({
                'gene': gene,
                'gene_row_index': i,
                'gene_y_coord': gene_y_coord,  # Where the gene actually is
                'label_y_coord': gene_y_coord  # Where we'll put the label (initially same)
            })
    
    # Apply repulsion to avoid overlapping labels
    if len(labeled_genes_data) > 1:
        total_height = len(actual_gene_order)
        
        # Sort by gene position
        labeled_genes_data.sort(key=lambda x: x['gene_y_coord'], reverse=True)
        
        # Adjust label positions to avoid overlap while keeping them in order
        for i in range(1, len(labeled_genes_data)):
            prev_label_y = labeled_genes_data[i-1]['label_y_coord']
            curr_label_y = labeled_genes_data[i]['label_y_coord']
            
            # If labels would overlap, move the current one down
            if prev_label_y - curr_label_y < min_label_spacing:
                labeled_genes_data[i]['label_y_coord'] = prev_label_y - min_label_spacing
        
        # Make sure no labels go below 0
        for i in range(len(labeled_genes_data)):
            if labeled_genes_data[i]['label_y_coord'] < 0:
                labeled_genes_data[i]['label_y_coord'] = 0
    
    # Add the labels and connecting lines

    for gene_data in labeled_genes_data:
        gene = gene_data['gene']
        gene_y = gene_data['gene_y_coord']  # Where the gene row actually is
        label_y = gene_data['label_y_coord']  # Where we'll place the label
        
        # Add the text label (potentially moved to avoid overlap)
        cm.ax_heatmap.text(
            x_outside,
            label_y,
            gene,
            fontsize=8,
            ha="left",
            va="center",
            transform=cm.ax_heatmap.transAxes,
            bbox=dict(facecolor='none', edgecolor='black', alpha=0.5, pad=1),
            clip_on=False
        )
        
        # Add a line connecting the label to the ACTUAL gene row
        cm.ax_heatmap.annotate(
            '',
            xy=(1.0, gene_y),  # Points to the actual gene row
            xytext=(x_outside - 0.002, label_y),  # Comes from the label position
            xycoords='axes fraction',
            textcoords='axes fraction',
            arrowprops=dict(arrowstyle='-', color='black', lw=0.5),
            clip_on=False
        )

        
def L_rank_genes_by_latent_time(
    adata,
    layer=None,
    key="latent_time",
    method="spearman",
    adjust_method="fdr_bh",
):
    """
    Rank genes by correlation with latent time, including adjusted p-values and filter flags.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    layer : str or None
        If specified, use this layer instead of `adata.X`.
    key : str
        Key in `adata.obs` to use as the latent time variable.
    method : str
        Correlation method: 'spearman' or 'pearson'.
    adjust_method : str
        Method for multiple testing correction (e.g., 'fdr_bh', 'bonferroni').

    Returns
    -------
    pd.DataFrame
        Columns: ['gene', 'correlation', 'pval', 'pval_adj', 'is_constant', 'is_all_nan']
    """
    import numpy as np
    import pandas as pd
    from scipy.stats import spearmanr, pearsonr
    from statsmodels.stats.multitest import multipletests

    if key not in adata.obs:
        raise ValueError(f"'{key}' not found in adata.obs.")

    latent = adata.obs[key].values
    X = adata.layers[layer] if layer and layer in adata.layers else adata.X
    X = X.toarray() if hasattr(X, "toarray") else X
    gene_names = adata.var_names

    corr_func = spearmanr if method == "spearman" else pearsonr

    corrs, pvals, is_const, is_nan = [], [], [], []

    for i in range(X.shape[1]):
        expr = X[:, i]
        constant = np.std(expr) == 0
        allnan = np.isnan(expr).all()

        is_const.append(constant)
        is_nan.append(allnan)

        if constant or allnan:
            corrs.append(np.nan)
            pvals.append(np.nan)
        else:
            corr, pval = corr_func(expr, latent)
            corrs.append(corr)
            pvals.append(pval)

    _, pvals_adj, _, _ = multipletests(pvals, method=adjust_method)

    result = pd.DataFrame({
        "gene": gene_names,
        "correlation": corrs,
        "pval": pvals,
        "pval_adj": pvals_adj,
        "is_constant": is_const,
        "is_all_nan": is_nan,
    })

    return result.sort_values("correlation", ascending=False).reset_index(drop=True)


def L_volcanoPlot_gene_latent_time_cor(
    cor_df,
    pval_key="pval_adj",
    r_thresh=0.3,
    p_thresh=0.05,
    figsize=(8, 6),
    title="Volcano Plot of Gene–Time Correlation",
    annotate_top=0,
):
    """
    Volcano-style plot of gene–latent time correlation with text repel.

    Parameters
    ----------
    cor_df : pd.DataFrame
        Must include 'gene', 'correlation', and the selected pval_key.
    pval_key : str
        Column name to use for p-values ('pval' or 'pval_adj').
    r_thresh : float
        Minimum absolute correlation for significance.
    p_thresh : float
        Maximum p-value threshold for significance.
    figsize : tuple
        Plot dimensions.
    title : str
        Plot title.
    annotate_top : int
        Number of lowest-p genes to label regardless of thresholds.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from adjustText import adjust_text
    
    df = cor_df.copy()
    df = df[df["correlation"].notna()].copy()
    if pval_key not in df.columns:
        raise ValueError(f"'{pval_key}' not found in DataFrame.")

    df["-log10(p)"] = -np.log10(df[pval_key])
    df["color"] = "gray"
    df.loc[(df["correlation"] > r_thresh) & (df[pval_key] < p_thresh), "color"] = "red"
    df.loc[(df["correlation"] < -r_thresh) & (df[pval_key] < p_thresh), "color"] = "blue"

    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(df["correlation"], df["-log10(p)"], c=df["color"], alpha=0.7, edgecolor="none")
    ax.axvline(x=r_thresh, color="black", linestyle="--", linewidth=1)
    ax.axvline(x=-r_thresh, color="black", linestyle="--", linewidth=1)
    ax.axhline(y=-np.log10(p_thresh), color="black", linestyle="--", linewidth=1)
    ax.set_xlabel("Correlation with latent time")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title(title)

    # Collect text labels for repel
    texts = []
    labelled_genes = set()

    # Label genes that pass both thresholds
    sig = df[(df[pval_key] < p_thresh) & (np.abs(df["correlation"]) > r_thresh)]
    for _, row in sig.iterrows():
        texts.append(ax.text(row["correlation"], row["-log10(p)"], row["gene"],
                             fontsize=8, ha="right" if row["correlation"] < 0 else "left"))
        labelled_genes.add(row["gene"])

    # Label top N by p-value, regardless of threshold, but skip genes already labelled
    if annotate_top > 0:
        top_genes = df.nsmallest(annotate_top, pval_key)
        for _, row in top_genes.iterrows():
            if row["gene"] not in labelled_genes:
                texts.append(ax.text(row["correlation"], row["-log10(p)"], row["gene"],
                                     fontsize=8, ha="right" if row["correlation"] < 0 else "left"))
                labelled_genes.add(row["gene"])


    # Apply repel adjustment
    if texts:
        adjust_text(texts, arrowprops=dict(arrowstyle="-", color="gray", lw=0.5),verbose=False)

    plt.tight_layout()
    plt.show()


