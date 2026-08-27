#!/usr/bin/env python3
import sys
import Common
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import ElasticNetCV
from sklearn.decomposition import PCA
from scipy import stats


def get_args(args):
    parser = argparse.ArgumentParser(description="Processing post-GWAS to identify gene-level effect",
                                    formatter_class=Common.CustomFormatter
        )
    parser.add_argument("-s", "--site", required=True,help="SNP-gene annotation file (must contain ID, Gene, R2, Type)")
    parser.add_argument("-g", "--gt", required=True, help="Genotype file (AA/AB/BB)")
    parser.add_argument("-p", "--phe", required=True,help="Phenotype file")
    parser.add_argument("-c", "--cov", required=False,help="Covariance file")
    parser.add_argument("-o","--out-prefix",type=str, default="Effect",help="Output prefix")
    parser.add_argument("-m", "--method", type=str, default="ElasticNet",help="Method of SNP-gene effect, ElasticNet or PCA")
    parser.add_argument("--ld-r2", type=float, default=0.9,help="LD r2 threshold for pruning ")
    parser.add_argument("--min-snps", type=int, default=1,help="Minimum SNPs per gene after pruning")
    parser.add_argument("--l1-ratio", type=float, default=0.5,help="ElasticNet l1_ratio")


    return parser.parse_args(args)


############################################
# Main
############################################
def main(args=None):
    args = get_args(args)

    # ---------- read input ----------
    map_df = pd.read_csv(args.site, sep="\t")
    if "R2" not in map_df.columns:
        raise ValueError("R2 column not found in site file")

    geno_raw = pd.read_csv(args.gt, sep="\t")
    sample_cols = geno_raw.columns[5:]

    for s in sample_cols:
        geno_raw[s] = geno_raw.apply(
            lambda r: gt_to_dosage(r[s], r["Ref"], r["Alt"]),
            axis=1
        )

    geno = geno_raw.set_index("ID")[sample_cols].T

    phe = pd.read_csv(args.phe, sep="\t", header=None)
    phe.columns = ["sample", "sample2", "trait"]

    geno,phe = filter_NA(geno,phe)

    gene_df = parpare_gene_scores(map_df,phe,geno,ld_r2=args.ld_r2,min_snps=args.min_snps,out_prefix=args.out_prefix,method=args.method,l1_ratio=args.l1_ratio)

    r2 = compute_partial_r2(phe, gene_df,args.cov)
    print(f"# gene effect result write to file: {args.out_prefix}.gene_effect.{args.method}.tsv...")
    r2.to_csv(f"{args.out_prefix}.gene_effect.{args.method}.tsv", sep="\t",index=False)
    #gene_df = cal_gene_scores(Xg,y,out_prefix=args.out_prefix,method=args.m,l1_ratio=args.l1_ratio)  


def filter_NA(geno, phe):
    # 1. Remove samples with missing phenotype
    phe = phe.dropna(subset=["trait"]).copy()

    # 2. Keep samples shared by genotype and phenotype
    common_samples = phe.loc[
        phe["sample"].isin(geno.index), "sample"
    ].tolist()

    # 3. Filter and reorder phenotype
    phe = (
        phe.set_index("sample")
           .loc[common_samples]
           .reset_index()
    )

    # 4. Filter and reorder genotype
    geno = geno.loc[common_samples]

    return geno, phe
    

def parpare_gene_scores(map_df,phe,geno,ld_r2,min_snps,out_prefix,method,l1_ratio):
    
    #print(map_df)
    phe2 = phe.set_index("sample")["trait"]

    # ---------- align samples ----------
    geno = geno.loc[phe2.index]
    y = phe2.values

    snp_r2 = map_df.set_index("ID")["R2"].to_dict()
    #snp_type = map_df.set_index("ID")["R2"].to_dict()

    gene_scores = {"ElasticNet": {}, "PCA": {}}
    gene_support = {}

    # ---------- gene loop ----------
    for gene, sub in map_df.groupby("Gene"):
        
        snps = [s for s in sub["ID"].unique() if s in geno.columns]

        #print(snps)
        if len(snps) < min_snps:
            print(f"#SNPs in {gene} is less than {min_snps}")
            continue

        Xg_df = geno[snps]

        #snps_ld = ld_prune(Xg_df, ld_r2)
        snps_ld = Xg_df.columns.tolist()

        if len(snps_ld) < min_snps:
            print(f"#LD-SNPs in {gene} is less than {min_snps}")
            continue

        Xg = Xg_df[snps_ld].values

        weights1 = np.array([np.sqrt(snp_r2.get(s)) for s in snps_ld])
        weights2 = np.array([get_snp_type_weight(sub, s) for s in snps_ld])

        weights = weights1 * weights2
        Xg = Xg * weights


        scaler = StandardScaler()
        Xg_std = scaler.fit_transform(Xg)
        # Elastic Net
        if method == "ElasticNet":
            G_enet, coef = elasticnet_gene_score(Xg_std, y, l1_ratio)
            if G_enet is not None:
                gene_scores["ElasticNet"][gene] = G_enet
                gene_support.setdefault(gene, []).append("ElasticNet")
            else:
                print(f"#{gene} can't cal ElasticNet")
                #print(Xg)
        # PCA
        if method == "PCA":
            G_pca, var_exp = pca_gene_score(Xg_std)
            gene_scores["PCA"][gene] = G_pca
            if var_exp > 0.2:
                gene_support.setdefault(gene, []).append("PCA")
    
    
    # ---------- output ----------
    def zscore_df(df):
        return pd.DataFrame(
            StandardScaler().fit_transform(df),
            index=df.index,
            columns=df.columns
        )

    df = zscore_df(pd.DataFrame(gene_scores[method], index=geno.index))
    print(f"# gene score result write to file: {out_prefix}.gene_score.{method}.tsv...")
    df.to_csv(f"{out_prefix}.gene_score.{method}.tsv", sep="\t")

    return(df)
        #return Xg,y
        

def get_snp_type_weight(gene_sub_df, snp_id):

    snp_types = gene_sub_df.loc[gene_sub_df['ID'] == snp_id, 'Type'].tolist()

    if len(snp_types) == 1 and "eQTL" in snp_types: # trans eQTL : 1
        return 1

    if "NoSyn" in snp_types:
        w = 2
    else:
        w = 1

    if "eQTL" in snp_types: # cis eQTL: 1.5
        w += 1.5

    return w


def compute_partial_r2(phe_df, gene_score_df, cov_file=None, standardize=True):
    """
    Compute marginal or partial R² for each gene.

    Parameters
    ----------
    phe_df : pd.DataFrame
        Phenotype DataFrame with columns ["sample", "trait"] or index as sample.
        Example:
            sample     trait
            3W1-1     -1705116.2
            3W1-10    -1714850.6
    gene_score_df : pd.DataFrame
        Gene score matrix (samples x genes), index = sample
    cov_file : str or None
        Covariate file path (tab-delimited), first column = sample, rest = PCs
        If None, compute marginal R²
    standardize : bool
        Whether to z-score gene scores and covariates

    Returns
    -------
    pd.DataFrame
        Columns: ["gene", "R2", "type"]
    """

    # --------- 处理 phenotype ---------
    phe = phe_df.copy()
    if "sample" in phe.columns and "trait" in phe.columns:
        phe = phe.set_index("sample")["trait"]
    else:
        # 如果已经是 series 或 index = sample
        phe = phe.iloc[:, 0] if phe.ndim > 1 else phe

    # --------- 读取 covariates ---------
    covar_df = None
    if cov_file is not None:
        covar_df = pd.read_csv(cov_file, sep="\t", header=None)
        n_pc = covar_df.shape[1] - 1
        covar_df.columns = ["sample"] + [f"PC{i+1}" for i in range(n_pc)]
        covar_df = covar_df.set_index("sample")

    # --------- 样本对齐 ---------
    samples = set(phe.index) & set(gene_score_df.index)
    if covar_df is not None:
        samples &= set(covar_df.index)
    samples = sorted(samples)

    if len(samples) == 0:
        raise ValueError("No common samples found among phenotype, gene_score, and covariates!")

    y = phe.loc[samples].values
    gene_score_df = gene_score_df.loc[samples]
    if covar_df is not None:
        covar_df = covar_df.loc[samples]

    # --------- baseline matrix ---------
    if covar_df is not None:
        X_base = covar_df.values
        if standardize:
            X_base = StandardScaler().fit_transform(X_base)
        r2_type = "partial"
    else:
        X_base = np.ones((len(y), 1))  # intercept only
        r2_type = "marginal"

    # --------- baseline RSS ---------
    base_model = LinearRegression()
    base_model.fit(X_base, y)
    y_hat0 = base_model.predict(X_base)
    rss0 = np.sum((y - y_hat0) ** 2)

    # --------- per-gene R² ---------
    results = []
    for gene in gene_score_df.columns:
        G = gene_score_df[[gene]].values
        if standardize:
            G = StandardScaler().fit_transform(G)

        X_full = np.hstack([X_base, G])
        model = LinearRegression()
        model.fit(X_full, y)
        y_hat1 = model.predict(X_full)
        rss1 = np.sum((y - y_hat1) ** 2)

        r2 = (rss0 - rss1) / rss0
        results.append((gene, r2, r2_type))

    return pd.DataFrame(results, columns=["gene", "R2", "type"]).sort_values(by="R2", ascending=False).reset_index(drop=True)



def compute_partial_r2_permutation(
    phe_df,
    gene_score_df,
    cov_file=None,
    standardize=True,
    n_perm=1000,
    random_state=0
):
    """
    Compute partial/marginal R² for each gene with permutation p-value.

    Parameters
    ----------
    phe_df : pd.DataFrame
        Phenotype DataFrame with columns ["sample", "trait"] or index as sample
    gene_score_df : pd.DataFrame
        Gene score matrix (samples x genes), index = sample
    cov_file : str or None
        Covariate file path (tab-delimited). First column = sample
    standardize : bool
        Whether to z-score gene scores and covariates
    n_perm : int
        Number of permutations for p-value
    random_state : int
        Random seed for reproducibility

    Returns
    -------
    pd.DataFrame
        Columns: ["gene", "R2", "type", "rank", "p_perm"], sorted by R2 descending
    """

    rng = np.random.default_rng(random_state)

    # --------- phenotype ---------
    phe = phe_df.copy()
    if "sample" in phe.columns and "trait" in phe.columns:
        phe = phe.set_index("sample")["trait"]
    else:
        phe = phe.iloc[:, 0] if phe.ndim > 1 else phe

    # --------- covariates ---------
    covar_df = None
    if cov_file is not None:
        covar_df = pd.read_csv(cov_file, sep="\t", header=None)
        n_pc = covar_df.shape[1] - 1
        covar_df.columns = ["sample"] + [f"PC{i+1}" for i in range(n_pc)]
        covar_df = covar_df.set_index("sample")

    # --------- 样本对齐 ---------
    samples = set(phe.index) & set(gene_score_df.index)
    if covar_df is not None:
        samples &= set(covar_df.index)
    samples = sorted(samples)

    if len(samples) == 0:
        raise ValueError("No common samples found among phenotype, gene_score, and covariates!")

    y = phe.loc[samples].values
    gene_score_df = gene_score_df.loc[samples]
    if covar_df is not None:
        covar_df = covar_df.loc[samples]

    # --------- baseline matrix ---------
    if covar_df is not None:
        X_base = covar_df.values
        if standardize:
            X_base = StandardScaler().fit_transform(X_base)
        r2_type = "partial"
    else:
        X_base = np.ones((len(y), 1))
        r2_type = "marginal"

    # --------- baseline RSS ---------
    base_model = LinearRegression()
    base_model.fit(X_base, y)
    y_hat0 = base_model.predict(X_base)
    rss0 = np.sum((y - y_hat0) ** 2)

    # --------- per-gene R² + permutation p ---------
    results = []
    for gene in gene_score_df.columns:
        G = gene_score_df[[gene]].values
        if standardize:
            G = StandardScaler().fit_transform(G)

        # R²
        X_full = np.hstack([X_base, G])
        model = LinearRegression()
        model.fit(X_full, y)
        y_hat1 = model.predict(X_full)
        rss1 = np.sum((y - y_hat1) ** 2)
        r2 = (rss0 - rss1) / rss0

        # permutation test
        perm_r2 = []
        for _ in range(n_perm):
            y_perm = rng.permutation(y)
            model.fit(X_full, y_perm)
            y_hat_perm = model.predict(X_full)
            rss1_perm = np.sum((y_perm - y_hat_perm) ** 2)
            perm_r2.append((rss0 - rss1_perm) / rss0)

        p_perm = (np.sum(np.array(perm_r2) >= r2) + 1) / (n_perm + 1)

        results.append((gene, r2, r2_type, p_perm))

    df = pd.DataFrame(results, columns=["gene", "R2", "type", "p_perm"])
    df = df.sort_values(by="R2", ascending=False).reset_index(drop=True)
    df["rank"] = np.arange(1, len(df) + 1)
    return df



def plot_cumulative_r2(
    y,
    gene_score_df,
    covar_df,
    gene_order,
    max_genes=30
):
    """
    gene_order: genes sorted by importance
    """

    y = np.asarray(y)
    X_base = covar_df.values

    base_model = LinearRegression().fit(X_base, y)
    R2_base = r2_score(y, base_model.predict(X_base))

    cum_r2 = []

    for k in range(1, max_genes + 1):
        genes_k = gene_order[:k]
        X = np.hstack([X_base, gene_score_df[genes_k].values])

        model = LinearRegression().fit(X, y)
        R2 = r2_score(y, model.predict(X))
        cum_r2.append(R2 - R2_base)

    plt.figure(figsize=(5, 4))
    plt.plot(range(1, max_genes + 1), cum_r2, marker="o")
    plt.xlabel("Number of genes")
    plt.ylabel("Cumulative explained variance (R²)")
    plt.title("Cumulative gene-level explained variance")
    plt.tight_layout()
    plt.show()




############################################
# Genotype conversion
############################################
def gt_to_dosage(gt, ref, alt):
    if gt == ref + ref:
        return 0
    elif gt in (ref + alt, alt + ref):
        return 1
    elif gt == alt + alt:
        return 2
    else:
        return np.nan


############################################
# LD pruning (gene-level)
############################################
def ld_prune(geno_df, r2_thresh=0.9):
    snps = geno_df.columns.tolist()
    if len(snps) <= 1:
        return snps

    corr = geno_df.corr().values
    r2 = corr ** 2

    kept = []
    for i, snp in enumerate(snps):
        keep = True
        for kept_snp in kept:
            j = snps.index(kept_snp)
            if r2[i, j] > r2_thresh:
                keep = False
                break
        if keep:
            kept.append(snp)

    return kept


############################################
# Elastic Net gene score
############################################
def elasticnet_gene_score(Xg_std, y, l1_ratio):
    model = ElasticNetCV(
        l1_ratio=l1_ratio,
        cv=5,
        max_iter=10000
    )
    model.fit(Xg_std, y)
    coef = model.coef_

    if np.all(coef == 0):
        return None, None

    Gg = Xg_std @ coef
    return Gg, coef


############################################
# PCA gene score
############################################
def pca_gene_score(Xg_std):
    pca = PCA(n_components=1)
    pc1 = pca.fit_transform(Xg_std)[:, 0]
    var_exp = pca.explained_variance_ratio_[0]
    return pc1, var_exp


############################################
# SKAT-like kernel test
############################################
def skat_gene_test(Xg, y):
    y = y - np.mean(y)
    K = Xg @ Xg.T
    Q = y.T @ K @ y

    trace_K = np.trace(K)
    trace_K2 = np.trace(K @ K)

    if trace_K2 == 0:
        return np.nan, np.nan

    df = 2 * trace_K**2 / trace_K2
    scale = trace_K2 / (2 * trace_K)
    pval = stats.chi2.sf(Q / scale, df)

    return Q, pval


if __name__ == "__main__":
    main()
