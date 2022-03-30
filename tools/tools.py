def pval_to_star(pval, GSEA_threshold=True):
    """
    This function returns star-formatted p-value.
    
    Parameters
    ----------
    pval : float
        P-value
        
    Returns
    ----------
    Star-formatted p-value.
    """
    if pval < 0.0001:
        return "****"
    elif pval < 0.001:
        return "***"
    elif pval < 0.01:
        return "**"
    elif pval < 0.05:
        return "*"
    elif pval < 0.25 and GSEA_threshold:
        return "."
    else:
        return " "

def QC_histogram(adata, var_names=["total_counts", "n_genes_by_counts", "pct_counts_mt"], split_by=None):    
    """
    This is a helper function that draws three histograms with distribution of
    main scRNA-Seq experiment QC metrics.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    var_names : list, optional
        List with variables of interest
        
    Returns
    ----------
    Figure with three histograms.
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    fig, axs = plt.subplots(ncols=len(var_names), figsize=(len(var_names) * 4, 3))
    
    for i in range(len(var_names)):
        if split_by == None:
            sns.histplot(adata.obs[var_names[i]], ax=axs[i])
        else:
            color_counter = 0
            for group in set(adata.obs[split_by]):
                sns.histplot(adata[adata.obs[split_by] == group].obs[var_names[i]], ax=axs[i],
                             label=group, color=sns.color_palette()[color_counter % 10])
                color_counter += 1
        if i == 0:
            axs[i].set_ylabel("Number of cells")
        else:
            axs[i].set_ylabel("")
        axs[i].set_xlabel(var_names[i])
        if split_by != None:
            axs[i].legend()
    fig.tight_layout(pad=1)
    return fig

    
def get_count_pseudobulk(adata, layer=None, split_by=None, use_raw=False):
    """
    Returns DataFrame with pseudo-bulks generated based on AnnData object.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with single cell expressions.
    layer : str, optional
        Name of the layer with count data.
    split_by : str, optional
        Name of column in `adata.obs` with variable to split.
    use_raw : bool, optional
        If it's needed to use `adata.raw`

    Returns
    ----------
    DataFrame with pseudo-bulk.
    
    """
    if split_by is None:
        if use_raw:
            df_res = pd.DataFrame({"counts" : np.matrix(adata.raw.X.sum(axis=0)).A[0]})
        elif layer is None:
            df_res = pd.DataFrame({"counts" : np.matrix(adata.X.sum(axis=0)).A[0]})
        else:
            df_res = pd.DataFrame({"counts" : np.matrix(adata.layers[layer].sum(axis=0)).A[0]})
    else:
        df_res = pd.DataFrame()
        for splitter in set(adata.obs[split_by]):
            if use_raw:
                df_res[splitter] = np.matrix(adata[adata.obs[split_by] == splitter].raw.X.sum(axis=0)).A[0]
            elif layer is None:
                df_res[splitter] = np.matrix(adata[adata.obs[split_by] == splitter].X.sum(axis=0)).A[0]
            else:
                df_res[splitter] = np.matrix(adata[adata.obs[split_by] == splitter].layers[layer].sum(axis=0)).A[0]
    if use_raw:
        df_res.index = adata.raw.var_names
    else:
        df_res.index = adata.var_names
    return df_res

def get_beautiful_cmap(initial_cmap="Reds", grey_intensity=0.2, color_intencity=0.1):
    """
    Returns color map for visualization of gene expression on UMAPs. Color
    map will starts from grey, not from white.
    
    Parameters
    ----------
    initial_cmap : str, optional
        What color map will be the base for novel color map.
    grey_intensity : float, optional
        What intensity of grey should be at the start of color map.
    color_intencity : float, optional
        What intensity of color should be after grey at color map

    Returns
    ----------
    ListedColormap object.
    
    """
    from matplotlib import cm
    from matplotlib.colors import ListedColormap
    import numpy as np
    
    cm_color = cm.get_cmap(initial_cmap, 128)
    cm_grey = cm.get_cmap("Greys", 128)
    
    c = ListedColormap(
        np.vstack(
            (cm_grey(np.linspace(0.2, 0.2, 1)),
             cm_color(np.linspace(0.1, 1, 128)))
    ))
    
    return c


def find_doublets(adata, layer=None, batch_key=None, score_key="doublet_score", remove_doublets=False,
                  status_key="predicted_doublet"):
    """
    This function finds doublets in AnnData object using Scrublet (https://github.com/swolock/scrublet).
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    layer : str, optional
        Name of layer with count data.
    batch_key : str, optional
        Name of column in adata.obs with an information about batch. By default it considered that
        data consists only of one batch.
    score_key : str, optional
        Name of the column in adata.obs doublet score would be placed to. By default it is "doublet_score"
    status_key : str, optional
        Name of the column in adata.obs prediction status would be placed to. By default it is
        "predicted_doublet"
    delete_doublets : bool, optional
        If found doublets should be removed. False by default.

    Returns
    ----------
    None
    """
    import scrublet as scr
    import pandas as pd
    
    doublet_df_index = []
    doublet_df = {
        "doublet_score" : [],
        "predicted_doublet" : []
    }
    if batch_key is None:
        batch_list = [True]
    else:
        batch_list = set(adata.obs[batch_key])
    for batch in batch_list:
        if batch_key is None:
            doublet_df_index += list(adata.obs.index)
            if layer is None:
                scrub = scr.Scrublet(adata.X)
            else:
                scrub = scr.Scrublet(adata.layers[layer])
        else:
            print(f"Batch {batch} in progress...")
            doublet_df_index += list(adata[adata.obs[batch_key] == batch].obs.index)
            if layer is None:
                scrub = scr.Scrublet(adata[adata.obs[batch_key] == batch].X)
            else:
                scrub = scr.Scrublet(adata[adata.obs[batch_key] == batch].layers[layer])
        doublet_scores, predicted_doublets = scrub.scrub_doublets()
        doublet_df["doublet_score"] += list(doublet_scores)
        doublet_df["predicted_doublet"] += list(predicted_doublets)
    doublet_df = pd.DataFrame(doublet_df)
    doublet_df.index = doublet_df_index
    doublet_df = doublet_df.loc[adata.obs.index]
    adata.obs[score_key] = doublet_df["doublet_score"]
    adata.obs[status_key] = doublet_df["predicted_doublet"]
    print(f"{sum(adata.obs[status_key])} doublets (of {len(adata)}) were found.")
    if remove_doublets:
        adata = adata[np.logical_not(adata.obs[status_key])]
    adata.obs[status_key] = doublet_df["predicted_doublet"].astype(str)
    adata.obs[status_key] = adata.obs[status_key].astype("category")
        
def get_cluster_proportions(adata,
                            cluster_key="cluster_final",
                            sample_key="replicate",
                            drop_values=None):
    """
    Input
    =====
    adata : AnnData object
    cluster_key : key of `adata.obs` storing cluster info
    sample_key : key of `adata.obs` storing sample/replicate info
    drop_values : list/iterable of possible values of `sample_key` that you don't want
    
    Returns
    =======
    pd.DataFrame with samples as the index and clusters as the columns and 0-100 floats
    as values
    """
    
    adata_tmp = adata.copy()
    sizes = adata_tmp.obs.groupby([cluster_key, sample_key]).size()
    props = sizes.groupby(level=1).apply(lambda x: 100 * x / x.sum()).reset_index() 
    props = props.pivot(columns=sample_key, index=cluster_key).T
    props.index = props.index.droplevel(0)
    props.fillna(0, inplace=True)
    
    if drop_values is not None:
        for drop_value in drop_values:
            props.drop(drop_value, axis=0, inplace=True)
    return props

def get_cluster_sizes(adata,
                      cluster_key="cluster_final",
                      sample_key="replicate",
                      drop_values=None):
    """
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    cluster_key : str, optional
        Key of `adata.obs` storing cluster info.
    sample_key : str, optional
        Key of `adata.obs` storing sample/replicate info.
    drop_values : list, optional
        Possible values of `sample_key` that you don't want
        
    Returns
    ----------
    pd.DataFrame with samples as the index and clusters as the columns
    """
    
    adata_tmp = adata.copy()
    sizes = adata_tmp.obs.groupby([cluster_key, sample_key]).size()
    sizes = sizes.groupby(level=1).apply(lambda x: x).reset_index()
    sizes = sizes.pivot(columns=sample_key, index=cluster_key).T
    sizes.index = sizes.index.droplevel(0)
    sizes.fillna(0, inplace=True)
    
    if drop_values is not None:
        for drop_value in drop_values:
            sizes.drop(drop_value, axis=0, inplace=True)
    return sizes

def run_conos_pagoda2(adata, batch_key, tmp_dir, script_path, layer=None, space="PCA", only_read_outputs=False,
                      n_hvg=3000, n_comps=30):
    """
    It's a conos wrapper.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    batch_key : str
        String with name of column with batch annotation.
    tmp_dir : str
        Path to empty tmp folder.
    layer : str, optional
        Name of layer with count data.
    space : str, optional (`PCA` or `CCA`)
        Alignment space.
    only_read_outputs : str, optional
        Technical variable.
    n_hvg : int, optional
        Number of highly variable genes.
    script_path : str, optional
        Path to a wrapper script.
        
    Returns
    ----------
    Aligned AnnData (Annotated data matrix) object.
    """
    import anndata as ad
    import pandas as pd
    import subprocess
    from scipy import io
    
    if tmp_dir[-1] != "/":
        tmp_dir += "/"
    adata_tmp = adata.copy()
    if not(layer is None):
        adata_tmp.X = adata_tmp.layers[layer].copy()
    adata_tmp = ad.AnnData(
        X=adata_tmp.X,
        obs=pd.DataFrame(adata_tmp.obs[batch_key]),
        var=pd.DataFrame(index=adata_tmp.var_names)
    )
    if not only_read_outputs:
        print("(1/3) Writing .h5ad-file")
        path_to_adata = tmp_dir + "adata.h5ad"
        adata_tmp.write_h5ad(path_to_adata)

        print("(2/3) Running wrapper (it may takes a long time)")
        subprocess.call(["Rscript", "--vanilla", script_path, tmp_dir, batch_key, space, str(n_hvg), str(n_comps)])
    
    print("(3/3) Reading conos outputs")
    umap_df = pd.read_csv(tmp_dir + "embedding.csv")
    pca_df = pd.read_csv(tmp_dir + "pca.csv")
    pseudopca_df = pd.read_csv(tmp_dir + "pseudo_pca.csv")
    graph_conn_mtx = io.mmread(tmp_dir + "graph_connectivities.mtx")
    graph_dist_mtx = io.mmread(tmp_dir + "graph_distances.mtx")
    counts_mtx = io.mmread(tmp_dir + "count_matrix.mtx")
    
    adata_tmp.layers["counts"] = adata_tmp.X.copy()
    adata_tmp.obsm["X_umap"] = umap_df.values
    adata_tmp.obsm["X_pca"] = pca_df.values
    adata_tmp.obsm["X_pseudo_pca"] = pseudopca_df.values
    adata_tmp.uns["neighbors"] = dict(connectivities=graph_conn_mtx.tocsr(), distances=graph_dist_mtx.tocsr())
    adata_tmp.uns["neighbors"]["params"] = dict(n_pcs=30, use_rep="X_pca", metric="cosine", method="umap", n_neighbors=15)
    adata_tmp.X = counts_mtx.tocsr()
    adata_tmp.var_names_make_unique()
    adata_tmp.obs = adata.obs.copy()
    adata_tmp.write_h5ad(tmp_dir + "adata.h5ad")
    
    return adata_tmp

def correct_all_pvalues(df, correction_method="bonferroni"):
    from statsmodels.stats.multitest import multipletests
    return pd.DataFrame(
        data=np.reshape(multipletests(np.ravel(target_df.values), method=correction_method)[1], np.shape(target_df.values)),
        index=target_df.index,
        columns=target_df.columns
    )

def get_DESeq_multigroup(adata, split_by="sample_id", grouping="condition"):
    """
    It's a DESeq2 wrapper for DE testing within three or more groups.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    split_by : str, optional
        String with name of column with batch annotation.
    grouping : str, optional
        Path to empty tmp folder.
        
    Returns
    ----------
    Dictionary with DESeq2 results.
    """
    import pandas as pd
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri, Formula
    pandas2ri.activate()
    from rpy2.robjects.packages import importr
    deseq = importr("DESeq2")
    BiocGenerics = importr("BiocGenerics")
    
    to_dataframe = robjects.r('function(x) data.frame(x)')
    
    class py_DESeq2:
        def __init__(self, count_matrix, design_matrix, design_formula, gene_column='id'):
            try:
                assert gene_column in count_matrix.columns, 'Wrong gene id column name'
                gene_id = count_matrix[gene_column]
            except AttributeError:
                sys.exit('Wrong Pandas dataframe?')

            self.dds = None
            self.deseq_result = None
            self.resLFC = None
            self.comparison = None
            self.normalized_count_matrix = None
            self.gene_column = gene_column
            self.gene_id = count_matrix[self.gene_column]
            self.count_matrix = robjects.conversion.py2rpy(count_matrix.drop(gene_column,axis=1))
            self.design_matrix = robjects.conversion.py2rpy(design_matrix)
            self.design_formula = Formula(design_formula)

        def run_deseq(self, **kwargs):
            self.dds = deseq.DESeqDataSetFromMatrix(
                countData=self.count_matrix, 
                colData=self.design_matrix,
                design=self.design_formula
            )
            self.dds = deseq.DESeq(self.dds, **kwargs)
            self.normalized_count_matrix = BiocGenerics.counts(self.dds, normalized=True)

        def get_deseq_result(self, **kwargs):
            self.comparison = deseq.resultsNames(self.dds)
            self.deseq_result = deseq.results(self.dds, **kwargs)
            self.deseq_result = to_dataframe(self.deseq_result)
            self.deseq_result = robjects.conversion.rpy2py(self.deseq_result)
            self.deseq_result[self.gene_column] = self.gene_id.values
    
    expression_matrix = get_count_pseudobulk(adata, layer="counts", split_by=split_by).astype("int32")
    mapping = dict(zip(
        adata.obs[[split_by, grouping]].drop_duplicates()[split_by],
        adata.obs[[split_by, grouping]].drop_duplicates()[grouping]
    ))
    conditions = pd.DataFrame({"treatment" : [mapping[sample] for sample in expression_matrix.columns]},
                              index=expression_matrix.columns)
    expression_matrix["id"] = expression_matrix.index
    expression_matrix.index = list(range(len(expression_matrix)))
    
    DESeq = py_DESeq2(expression_matrix, conditions, design_formula="~ treatment")
    DESeq.run_deseq()
    
    results = {
        "DE" : {},
        "counts" : get_count_pseudobulk(adata, layer="counts", split_by=split_by).astype("int32"),
        "normalized counts" : pd.DataFrame(
            data=DESeq.normalized_count_matrix,
            index=expression_matrix["id"],
            columns=expression_matrix.columns[:-1]
        ),
        "mapping" : mapping
    }
    
    for first_cond in set(adata.obs[grouping]):
        results["DE"][first_cond] = {}
        for second_cond in set(adata.obs[grouping]):
            if first_cond != second_cond:
                DESeq.get_deseq_result(contrast=robjects.vectors.StrVector(["treatment", first_cond, second_cond]))
                results["DE"][first_cond][second_cond] = DESeq.deseq_result
    
    return results

def get_DESeq(adata, group, label="cell_type_l2", split_by="sample_id", reference=False):
    """
    It's a DESeq2 wrapper for DE testing between two groups.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    group : str
        Name of target group.
    label : str, optional
        Name of column in `adata.obs` with target group labeling.
    split_by : str, optional
        Name of column in `adata.obs` with sample labeling.
    reference : str, optional
        Name of column in `adata.obs` with reference group labeling.
        If False use all non-target groups as reference.
        
    Returns
    ----------
    Dictionary with DESeq2 results.
    """
    import pandas as pd
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri, Formula
    pandas2ri.activate()
    from rpy2.robjects.packages import importr
    deseq = importr("DESeq2")
    BiocGenerics = importr("BiocGenerics")
    
    to_dataframe = robjects.r('function(x) data.frame(x)')
    
    class py_DESeq2:
        def __init__(self, count_matrix, design_matrix, design_formula, gene_column='id'):
            try:
                assert gene_column in count_matrix.columns, 'Wrong gene id column name'
                gene_id = count_matrix[gene_column]
            except AttributeError:
                sys.exit('Wrong Pandas dataframe?')

            self.dds = None
            self.deseq_result = None
            self.resLFC = None
            self.comparison = None
            self.normalized_count_matrix = None
            self.gene_column = gene_column
            self.gene_id = count_matrix[self.gene_column]
            self.count_matrix = robjects.conversion.py2rpy(count_matrix.drop(gene_column,axis=1))
            self.design_matrix = robjects.conversion.py2rpy(design_matrix)
            self.design_formula = Formula(design_formula)

        def run_deseq(self, **kwargs):
            self.dds = deseq.DESeqDataSetFromMatrix(
                countData=self.count_matrix, 
                colData=self.design_matrix,
                design=self.design_formula
            )
            self.dds = deseq.DESeq(self.dds, **kwargs)
            self.normalized_count_matrix = BiocGenerics.counts(self.dds, normalized=True)

        def get_deseq_result(self, **kwargs):
            self.comparison = deseq.resultsNames(self.dds)
            self.deseq_result = deseq.results(self.dds, **kwargs)
            self.deseq_result = to_dataframe(self.deseq_result)
            self.deseq_result = robjects.conversion.rpy2py(self.deseq_result)
            self.deseq_result[self.gene_column] = self.gene_id.values
    
    expr = get_count_pseudobulk(adata[adata.obs[label] == group], layer="counts", split_by=split_by).astype("int32")
    if not reference:
        expr_ref = get_count_pseudobulk(adata[adata.obs[label] != group], layer="counts", split_by=split_by).astype("int32")
    else:
        expr_ref = get_count_pseudobulk(adata[adata.obs[label] == reference], layer="counts", split_by=split_by).astype("int32")
    expr_ref.columns = expr_ref.columns + "_control"
    expression_matrix = expr.join(expr_ref)
    conditions = pd.DataFrame({"treatment" : ["tagret"] * len(expr.columns) + ["control"] * len(expr_ref.columns)}, index=expression_matrix.columns)
    expression_matrix["id"] = expression_matrix.index
    expression_matrix.index = list(range(len(expression_matrix)))
    
    DESeq = py_DESeq2(expression_matrix, conditions, design_formula="~ treatment")
    DESeq.run_deseq()
    DESeq.get_deseq_result()
    
    return DESeq.deseq_result

def get_GSEA(DESeq_results, genesets, outdir):
    """
    Performs GSEA analysis based on DESeq results.
    
    Parameters
    ----------
    DESeq_results : DataFrame
        Output of get_DESeq() function.
    genesets : dict
        Dictionary with gene sets for GSEA analysis.
    outdir : str
        Path to directory with outputs.
        
    Returns
    ----------
    DataFrame with GSEA results.
    """
    import gseapy as gp
    
    DESeq_df = DESeq_results.dropna().sort_values(by="log2FoldChange", ascending=False)
    rnk = pd.DataFrame({0 : DESeq_df.id, 1 : DESeq_df.log2FoldChange})

    pre_res = gp.prerank(
        rnk=rnk,
        gene_sets=genesets,
        processes=35,
        permutation_num=100,
        outdir=outdir,
        format="pdf",
        seed=0,
        no_plot=True,
        verbose=True
    )

    results = pre_res.res2d.copy()
    return results

def compare_pcts(group1_pcts, group2_pcts):
    """
    BetaReg wrapper for comparison of percentage data.
    
    Parameters
    ----------
    group1_pcts : list
        List with fractions within group 1.
    group2_pcts : list
        List with fractions within group 2.
        
    Returns
    ----------
    p-value of regression coefficient.
    """
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri, Formula
    from rpy2.robjects.packages import importr
    import pandas as pd

    pandas2ri.activate()
    betareg = importr("betareg")
    base = importr("base")
    
    df = robjects.conversion.py2rpy(pd.DataFrame({
        "group" : ["A"] * len(group1_pcts) + ["B"] * len(group2_pcts),
        "pcts" : group1_pcts + group2_pcts
    }))
    regression = betareg.betareg(Formula("pcts ~ group"), data = df)

    return base.summary(regression)[0][0].flatten()[7]

def multiple_compare_pcts(comparisons, correction_method="bonferroni", pseudo_count=1e-5):
    """
    Performs multiple beta regressions with multiple testing correction.
    
    Parameters
    ----------
    comparisons : list
        List of tuples (group1_pcts, group2_pcts), see compare_pcts()
        description.
    correction_method : str, optional
        Name of multiple testing correction method (see
        statsmodels.stats.multitest.multipletests() description).
    pseudo_count : float, optional
        Pseudocount for dealing with zero values.
        
    Returns
    ----------
    p-values of comparisons.
    """
    from statsmodels.stats.multitest import multipletests
    ps = []
    for comparison in comparisons:
        group1_pcts = np.array(comparison[0]) + pseudo_count
        group2_pcts = np.array(comparison[1]) + pseudo_count
        if ((group1_pcts != pseudo_count).sum() + (group2_pcts != pseudo_count).sum()) == 0:
            ps.append(1)
        else:
            ps.append(compare_pcts(list(group1_pcts), list(group2_pcts)))
    if correction_method is None:
        return ps
    else:
        return multipletests(pvals=ps, method=correction_method)[1]

def multuple_celltype_condition_comparisons(df, groups_to_compare, correction_method="bonferroni",
                                            index=None, pseudo_count=1e-5):
    """
    Performs multiple beta regressions with multiple testing correction.
    
    Parameters
    ----------
    df : DataFrame
        DataFrame with cell fractions per sample. Columns are cell types,
        rows are samples named by sample name.
    groups_to_compare : list
        List of tuples ([samples_in_group_1], [samples_in_group_2]). It is needed
        to create `comparisons` argument for multiple_compare_pcts() function.
    correction_method : str, optional
        Name of multiple testing correction method (see
        statsmodels.stats.multitest.multipletests() description).
    index : list, optional
        List with names of comparisons.
    pseudo_count : float, optional
        Pseudocount for dealing with zero values.
        
    Returns
    ----------
    DataFrame with corrected p-values of comparisons.
    """
    import pandas as pd
    
    comparisons = []
    for celltype in df.columns:
        for comparison in groups_to_compare:
            comparisons.append(tuple([list(df[celltype][comparison[0]]), list(df[celltype][comparison[1]])]))
    ps_corrected = multiple_compare_pcts(comparisons, correction_method=correction_method, pseudo_count=pseudo_count)
    res = {}
    counter = 0
    for celltype in df.columns:
        res[celltype] = []
        for comparison in groups_to_compare:
            res[celltype].append(ps_corrected[counter])
            counter += 1
    return pd.DataFrame(res, index=index)