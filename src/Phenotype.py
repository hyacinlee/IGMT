#!/usr/bin/env python3 
import sys
import Common
import random
import time
import numpy as np
import pandas as pd
import seaborn as sns
import scipy.cluster.hierarchy as sch
from scipy.stats import ttest_ind
import argparse
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import BayesianRidge
from sklearn.impute import SimpleImputer
from joblib import Parallel, delayed



def get_args(args):
    parser = argparse.ArgumentParser(description="Processing phenotypic data, including clustering, dimensionality reduction, etc.",
                                    formatter_class=Common.CustomFormatter
        )
    parser.add_argument("-p","--phe", required=True, help="Path of the input phenotype file")
    parser.add_argument("-o","--out",help="Prefix of outfiles",type=str,default="out") 
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument("-g","--globals",help="Global PCA and automatically cluster dimensionality reduction",action='store_true')
    group.add_argument("-l","--locals",help="Local PCA and automatically cluster dimensionality reduction with specified traits")
    parser.add_argument("-i","--impute",required=False,choices=['iterative_bayes', 'iterative_rf'],help="impute method,choose from ['iterative_bayes', 'iterative_rf']")
    parser.add_argument("-c","--corr_threshold",help="Threshold of the cor bewteen two traits to divide into a same group,default=0.6",default=0.6,type=float)
    parser.add_argument("-z","--pc_threshold",help="Threshold of the PC1 bewteen two traits to divide into a same group,default=0.8",default=0.8,type=float)
    parser.add_argument("-m","--miss_threshold",help="The maximum allowable proportion of missing values. Those exceeding this proportion will be discarded,default=0.2",type=float,default=None)
    parser.add_argument("-f","--fig",help="Also output the groups corr-heatmap",action='store_true')
    parser.add_argument("-b","--bootstrap",help="Bootstrap times to run Cluster",default=500,type=int)
    parser.add_argument("--index",help="Add or change Index Name in output file",default="Sample",type=str)
    parser.add_argument("--mean",help="Calculate the average value of multiple columns using the specified grouping： name | groups",default=None,type=str)
    parser.add_argument("--group_ttest",help="Conduct a t-test using the specified group: name | groups",default=None,type=str)
    parser.add_argument("--stat",help="Commonly used metrics for line-based statistics",default=None,action='store_true')
    parser.add_argument("--cor",help="Cal and output cor bewteen two traits",default=None,action='store_true')
    parser.add_argument("--combine",help="Combine Column names",default=None)
    parser.add_argument("--combine_files",help="Files combined to -phe",default=None,nargs="*")
    parser.add_argument("--combine_method",help="Files combined method,choose from Base, Union or Intersection",default="Base")
    parser.add_argument("--subColumns",help="Extract the Columns in the file, in columns.",default=None,type=str)
    parser.add_argument("--subLines",help="Extract the Lines in the file, in lines.",default=None,type=str)
    parser.add_argument("--transpose",help="transpose line and columns and output",action='store_true')

    parser.add_argument("--top_k",help="Impute Optional: Maximum number of most correlated traits selected as predictors",default=10,type=int)
    parser.add_argument("--n_iter",help="Impute Optional: Maximum number of iterative imputation rounds ",default=5,type=int)
    parser.add_argument("--threads",help="Impute Optional: Number of parallel threads",default=24,type=int)
    parser.add_argument("--database", type=str,help="Impute Optional: Database matrix")

 
    parsed_args = parser.parse_args(args)
    
    return parsed_args
    

def main(args=None):

    args=get_args(args)

    phe_data=load_data(args.phe)

    if args.miss_threshold:
        phe_data = missing_values(phe_data, args.miss_threshold,args.out)


    if args.globals or args.locals:
        run_phe(phe_data,args)

    elif args.impute:
        database = None
        if args.database:
            database = pd.read_csv(args.database, sep="\t", index_col=0)
        impute_data(phe_data,output=args.out,database=database,method=args.impute,top_k=args.top_k,cor_threshold=args.corr_threshold,n_iter=args.n_iter,threads=args.threads)


    elif args.mean:
        run_matrix_mean(phe_data,args)

    elif args.group_ttest:
        run_matrix_group_ttest(phe_data,args)

    elif args.subColumns or args.subLines:
        run_matrix_sub(phe_data,args)

    elif args.transpose:
        run_matrix_transpose(phe_data,args)

    elif args.stat:
        row_statistics(phe_data,args.out)

    elif args.combine:
        lis=[args.phe]+args.combine_files
        run_combine(file_list=lis, on_col=args.combine, output=f"{args.out}.combine.tsv", method=args.combine_method)

    elif args.cor:
        cor_fig(phe_data,args.out)




def run_combine(file_list,on_col,output,method="Union",fillna="None", sep="\t",):

    if method not in ["Base", "Union", "Intersection"]:
        raise ValueError("method must be one of: Base, Union, Intersection")

    merged = None
    key_sets = []

    dfs = []
    for f in file_list:
        df = pd.read_csv(f, sep=sep, dtype=str)
        if on_col not in df.columns:
            raise ValueError(f"{f} there is no colum: {on_col}")

        dfs.append((f, df))
        key_sets.append(set(df[on_col].astype(str)))


    if method == "Base":
        keep_keys = key_sets[0]  
    elif method == "Union":
        keep_keys = set.union(*key_sets)
    elif method == "Intersection":
        keep_keys = set.intersection(*key_sets)
    else:
        raise ValueError("Invalid merge method")


    filtered_dfs = []
    for f, df in dfs:
        df = df[df[on_col].astype(str).isin(keep_keys)].copy()
        filtered_dfs.append(df)


    merged = filtered_dfs[0]
    for df in filtered_dfs[1:]:
        merged = pd.merge(merged, df, on=on_col, how="outer")


    merged = merged.fillna(fillna)

    if output:
        merged.to_csv(output, sep=sep, index=False)
    else:
        print(merged)

    return merged



def row_statistics(df,output):
    df_num = df.select_dtypes(include=[np.number])

    stats_df = pd.DataFrame({
        'sum': df_num.sum(axis=1),
        'mean': df_num.mean(axis=1),
        'min': df_num.min(axis=1),
        'max': df_num.max(axis=1),
        'median': df_num.median(axis=1),
        'std': df_num.std(axis=1),
        'count': df_num.count(axis=1),
        'nonzero': (df_num != 0).sum(axis=1)
    })
    stats_df = stats_df.reset_index()
    stats_df.rename(columns={'index': 'Info'}, inplace=True)
    stats_df.to_csv(f"{output}.stat_row.tsv",sep="\t",index=False)



def run_matrix_group_ttest(df,args):
    print("Conduct a t-test using the specified group")
    sample_to_group=Common.read_file(args.group_ttest,mode="dict",vals=[1],keys=[0])
    valid_samples = [s for s in df.columns if s in sample_to_group]
    df = df[valid_samples]

    groups = list(set(sample_to_group[s] for s in valid_samples))
    if len(groups) != 2:
        raise ValueError(f"Exactly two sets are needed, but there are now {len(groups)} groups: {groups}")

    g1, g2 = groups
    g1_samples = [s for s in valid_samples if sample_to_group[s] == g1]
    g2_samples = [s for s in valid_samples if sample_to_group[s] == g2]

    results = []
    for gene, row in df.iterrows():
        vals1 = row[g1_samples].astype(float).values
        vals2 = row[g2_samples].astype(float).values
        
        # T-test
        stat, pval = ttest_ind(vals1, vals2, equal_var=False, nan_policy="omit")
        
        g1_mean, g1_std = vals1.mean(), vals1.std(ddof=1)
        g2_mean, g2_std = vals2.mean(), vals2.std(ddof=1)
        
        # log2 fold change
        #log2fc = np.log2((g2_mean + 1) / (g1_mean + 1))
        
        results.append({
            "gene_id": gene,
            "pvalue": pval,
            #"log2FC": log2fc,
            f"{g1}_mean": g1_mean,
            f"{g1}_std": g1_std,
            f"{g2}_mean": g2_mean,
            f"{g2}_std": g2_std     
        })

    res_df = pd.DataFrame(results).sort_values("pvalue")
    res_df.to_csv(f"{args.out}.grouped.t-test.tsv",sep="\t",index=False)
    print(res_df)
    #return res_df



def run_matrix_mean(df,args):
    print("Calculate the average value of multiple columns using the specified grouping")
    sample_to_group=Common.read_file(args.mean,mode="dict",vals=[1],keys=[0])
    valid_samples = [s for s in df.columns if s in sample_to_group]
    df = df[valid_samples]
    df_renamed = df.rename(columns=sample_to_group)
    df_grouped = df_renamed.groupby(level=0, axis=1).mean()
    df_grouped.to_csv(f"{args.out}.grouped.tsv",sep="\t")
    #print(df_grouped)


def run_matrix_sub(phe_data,args):
    print("Get Sub df with samples or traits")
    df_sub=sub_sample_traits(phe_data,args.subColumns,args.subLines)
    df_sub.to_csv(f"{args.out}.sub_set.tsv",sep="\t")
    exit(0)

def run_matrix_transpose(phe_data,args):
    df = phe_data.T
    df.index.name = args.index
    df.to_csv(f"{args.out}.transpose.tsv",sep="\t")
    exit(0)


def run_phe(phe_data,args):

    #df_mis = missing_values(phe_data, args.miss_threshold,args.out)

    df_mis_inpute = phe_data

    if args.globals:
        best_result, change_history = bootstrap_group_traits(df_mis_inpute, args.corr_threshold, args.pc_threshold, args.bootstrap)
        df_change = pd.DataFrame(change_history, columns=['times', 'length', 'n50', 'uncluster', 'top20'])
        df_change.to_csv(f"{args.out}.best_history.tsv", sep='\t', index=False)
        
        groups, pc_values, group_dfs, group_pc_dfs = best_result
        output_cluster(df_mis_inpute, groups, pc_values, group_dfs, group_pc_dfs, args.out, args.fig)


    if args.locals:
        start_list = pd.read_csv(args.locals, sep='\t').iloc[:, 0].tolist()
        best_result, change_history = bootstrap_group_traits(df_mis_inpute, args.corr_threshold, args.pc_threshold, args.bootstrap, start_list=start_list)
        df_change = pd.DataFrame(change_history, columns=['times', 'length', 'n50', 'uncluster', 'top20'])
        df_change.to_csv(f"{args.out}.locals_best_history.tsv", sep='\t', index=False)

        groups, pc_values, group_dfs, group_pc_dfs = best_result
        output_cluster(df_mis_inpute, groups, pc_values, group_dfs, group_pc_dfs, args.out, args.fig, plot_bar=False)



    exit(0)

    '''
    if args.impute:
        df_mis_inpute=imputer_data(df_mis,output=args.out,method=args.impute,max_iter=args.max_iter, random_state=args.random_state,n_estimators=args.n_estimators,n_neighbors=args.n_neighbors)
        #print(df_mis_inpute)
        #imputer_data(df,method="iterative_rf",max_iter=5, random_state=0,n_estimators=100,n_neighbors=5)

    if args.globals:

        best_group_length = float('inf')
        best_n50 = -1
        best_uncluster = float('inf')
        best_top20 = -1
        result=[]
        change=[]

        for i in range(args.bootstrap):
            groups,pc_values,group_dfs,group_pc_dfs = group_traits(df_mis_inpute,args.corr_threshold,args.pc_threshold)

            group_length=len(groups)

            n50,top_20_ratio,uncluster = n50_from_2d_array(groups)

            is_better = False
    
            if group_length < best_group_length:
                is_better = True
            elif group_length == best_group_length:
                if n50 > best_n50:
                    is_better = True
                elif n50 == best_n50:
                    if uncluster < best_uncluster:
                        is_better = True
                    elif uncluster == best_uncluster:
                        if top_20_ratio > best_top20:
                            is_better = True

            if is_better:
                best_group_length = group_length
                best_n50 = n50
                best_uncluster = uncluster
                best_top20 = top_20_ratio
                result=[groups,pc_values,group_dfs,group_pc_dfs]

            change.append([i+1,best_group_length,best_n50,best_uncluster,best_top20])
            print(f"# Bootstrap {i+1} times，group lengths = {len(groups)} ;N50 = {n50} ; top_20_ratio = {top_20_ratio} ;  uncluster = {uncluster}")

        groups,pc_values,group_dfs,group_pc_dfs = result
        output(df_mis_inpute,groups,pc_values,group_dfs,group_pc_dfs,args.out,args.fig)
        df_change = pd.DataFrame(change, columns=['times', 'length', 'n50', 'uncluster', 'top20'])
        df_change.to_csv(f"{args.out}.best_history.tsv", sep='\t', index=False)


    if args.locals: #locals
        start_list=pd.read_csv(args.locals, sep='\t').iloc[:, 0].tolist()
        #for i in range(args.bootstrap):
        groups,pc_values,group_dfs,group_pc_dfs = group_traits(df_mis_inpute,args.corr_threshold,args.pc_threshold,start_list)
        output(df_mis_inpute,groups,pc_values,group_dfs,group_pc_dfs,args.out,args.fig)
    '''




def bootstrap_group_traits(df_mis_inpute, corr_threshold, pc_threshold, bootstrap_times, start_list=None):
    """
    Perform bootstrap grouping and select the best result.
    If start_list is provided, it will be shuffled each iteration.

    Returns:
        best_result: [groups, pc_values, group_dfs, group_pc_dfs]
        history: list of [iteration, best_group_length, best_n50, best_uncluster, best_top20]
    """
    best_group_length = float('inf')
    best_n50 = -1
    best_uncluster = float('inf')
    best_top20 = -1
    best_result = []
    history = []

    for i in range(bootstrap_times):
        if start_list is not None:
            # 打乱 start_list
            shuffled_start = start_list.copy()
            random.shuffle(shuffled_start)
            groups, pc_values, group_dfs, group_pc_dfs = group_traits(
                df_mis_inpute, corr_threshold, pc_threshold, shuffled_start
            )
        else:
            groups, pc_values, group_dfs, group_pc_dfs = group_traits(
                df_mis_inpute, corr_threshold, pc_threshold
            )

        group_length = len(groups)
        n50, top_20_ratio, uncluster = n50_from_2d_array(groups)

        # 判断是否比当前最佳更好
        is_better = (
            (group_length < best_group_length) or
            (group_length == best_group_length and (
                (n50 > best_n50) or
                (n50 == best_n50 and (
                    (uncluster < best_uncluster) or
                    (uncluster == best_uncluster and top_20_ratio > best_top20)
                ))
            ))
        )

        if is_better:
            best_group_length = group_length
            best_n50 = n50
            best_uncluster = uncluster
            best_top20 = top_20_ratio
            best_result = [groups, pc_values, group_dfs, group_pc_dfs]

        history.append([i+1, best_group_length, best_n50, best_uncluster, best_top20])

        print(
            f"# Bootstrap {i+1} times，group lengths = {group_length} ; "
            f"N50 = {n50} ; top_20_ratio = {top_20_ratio} ; uncluster = {uncluster}"
        )

    return best_result, history





def n50_from_2d_array(arr2d):
    lengths = [len(row) for row in arr2d]
    
    lengths_sorted = sorted(lengths, reverse=True)
    #print(lengths_sorted)
    top_20_percent_count = max(1, int(len(lengths_sorted) * 0.2))  # 至少 1 个
    top_20_sum = sum(lengths_sorted[:top_20_percent_count])
    total_sum = sum(lengths_sorted)
    top_20_ratio = top_20_sum / total_sum

    count_len1 = sum(1 for l in lengths if l == 1)

    total = sum(lengths_sorted)
    cum_sum = 0
    
    for l in lengths_sorted:
        cum_sum += l
        if cum_sum >= total / 2:
            n50 = l
            break
    
    return [n50,top_20_ratio,count_len1]

def sub_sample_traits(df,subColumns,subLines):
    #print(df.columns)
    #print(df)
    df.index = df.index.astype(str)
    if subColumns:
        trait_list=Common.read_file(subColumns,mode="list",vals=[0])
        print(f"# Keep {len(trait_list)} Columns.")
        missing_traits = [t for t in trait_list if t not in df.columns]
        sub_traits = [t for t in trait_list if t in df.columns]
        if missing_traits:
            print(f"The following traits were not found in the DataFrame and will be ignored: {missing_traits}")
        df = df.loc[:,sub_traits]
    if subLines:
        samples_list=Common.read_file(subLines,mode="list",vals=[0])
        print(f"# Keep {len(samples_list)} Lines.")
        df.index = df.index.astype(str)
        valid_samples = [str(s) for s in samples_list if s in df.index]
        df = df.loc[valid_samples]

    return df



def group_traits(df, corr_threshold=0.8, pc_threshold=0.8, start_list=None):
    corr_matrix = df.corr()  # 保留正负号

    all_traits = list(df.columns)
    random.shuffle(all_traits)
    remaining_traits = all_traits.copy()

    grouped_traits = []
    group_pc_values = []
    group_dfs = []
    group_pc_dfs = []

    start_queue = list(start_list) if start_list else []

    while remaining_traits:
        if start_list:   # local
            if start_queue:
                base_trait = start_queue.pop()
                if base_trait not in remaining_traits:
                    continue
                #print (f"#Cluster starts with {base_trait}")
            else:
                break
        else: # Global
            base_trait = remaining_traits.pop()

        group = [base_trait]

        if base_trait in remaining_traits:
            remaining_traits.remove(base_trait)

        #只考虑正相关性
        candidates = {
            trait for trait in remaining_traits
            if corr_matrix.loc[base_trait, trait] > corr_threshold
        }

        # 按 group 中已有性状的平均正相关性降序排序
        candidates = sorted(
            candidates,
            key=lambda trait: np.mean([corr_matrix.loc[trait, t] for t in group]),
            reverse=True
        )

        for trait in candidates:
            temp_group = group + [trait]
            temp_df = df[temp_group].dropna()

            if temp_df.shape[0] > 1:
                pc1_ratio, _ = compute_pca(temp_df)
                if pc1_ratio >= pc_threshold:
                    group.append(trait)
                    remaining_traits.remove(trait)

        group_df = df[group].dropna()
        if group_df.shape[1] == 0:
            continue

        pc1_ratio, pc1_df = compute_pca(group_df)
        grouped_traits.append(group)
        group_pc_values.append(pc1_ratio)
        group_dfs.append(group_df)
        group_pc_dfs.append(pc1_df)

    sorted_groups = sorted(
        zip(grouped_traits, group_pc_values, group_dfs, group_pc_dfs),
        key=lambda x: len(x[0]),
        reverse=True
    )

    return zip(*sorted_groups)



def output_cluster(df, groups, pc_values, group_dfs, group_pc_dfs, outprefix, fig, plot_bar=True):

    out_path = f"{outprefix}.Combine.Triat.txt"
    out = open(out_path, "w")

    df1 = pd.DataFrame(index=df.index)
    for i, pc1_scores in enumerate(group_pc_dfs, start=1):
        df1[f"Group{i}"] = pc1_scores["PC1"]

    for i, (group, pc1, dfg, dfpc) in enumerate(
            zip(groups, pc_values, group_dfs, group_pc_dfs), start=1):

        if len(group) > 2 and fig:
            cor_fig(dfg, f"{outprefix}.Combine.Group{i}.png")

        gg = ",".join(group)
        out.write(f"Group{i}\t{len(group)}\t{pc1:.2%}\t{group[0]}\t{gg}\n" )

    out.close()

    selected_traits = [trait for group in groups if len(group) >= 2 for trait in group]
    df3 = df[selected_traits]

    selected_traits_df4 = [group[0] for group in groups if len(group) == 1]
    df4 = df[selected_traits_df4]

    df1.to_csv(f"{outprefix}.PC1_scores.txt", sep="\t", na_rep="NA")
    df3.to_csv(f"{outprefix}.Trait.inGroup.txt", sep="\t", na_rep="NA")
    df4.to_csv(f"{outprefix}.Trait.unGroup.txt", sep="\t", na_rep="NA")

    if plot_bar:
        plot_group_bar(
            f"{outprefix}.Combine.Triat.txt",
            f"{outprefix}.Combine.Triat.pdf"
        )


def load_data(file_path):   
    df = pd.read_csv(file_path,sep="\t", index_col=0,dtype={0: str})  # 以第一列（样本名）作为索引
    df = df.loc[:, ~df.columns.str.startswith("Unnamed")]
    #df = pd.read_csv(file_path,sep="\t", index_col=0)
    #print(df)
    df.replace([np.inf, -np.inf], np.nan, inplace=True)  # 处理 inf

    return df


def missing_values(df,miss_threshold,out):
    print(f"The shape of input data: {df.shape}")
    df = df.select_dtypes(include=[np.number]) 
    missing_rates = df.isna().mean()
    columns_to_drop = missing_rates[missing_rates > miss_threshold].index
    df.drop(columns=columns_to_drop, inplace=True) 
    print(f"Deletes {len(list(columns_to_drop))} traits column with a deletion rate greater than {miss_threshold*100:.0f}% : {list(columns_to_drop)}")
    df.to_csv(f"{out}.dropmis.tsv",sep="\t", na_rep="NA")
    return df


#----------  ---------- ---------- impute ---------- ---------- ----------
def impute_data(data,output,database=None,method="iterative_rf",top_k=15,cor_threshold=0.4,n_iter=3,threads=24):
    """
    Impute missing values in a trait matrix.

    Parameters
    ----------
    data : pandas.DataFrame
        Target matrix (samples × traits).
    output : str
        Output file.
    database : pandas.DataFrame, optional
        Additional reference matrix. Samples should match the target matrix.
    method : {"iterative_rf", "iterative_bayes"}
        Imputation method.
    top_k : int
        Maximum number of predictor traits.
    cor_threshold : float
        Minimum absolute correlation for predictor selection.
    n_iter : int
        Number of iterative imputation rounds.
    threads : int
        Number of parallel threads.
    """

    start_time = time.time()

    df = data.copy()
    train_df = df.copy()

    if database is not None:
        common_samples = df.index.intersection(database.index)
        database = database.loc[common_samples]
        train_df = pd.concat([train_df, database], axis=1)

    print("# df shape:", df.shape)
    print("# train_df shape:", train_df.shape)

    print("# computing correlation matrix")
    train_fill = train_df.fillna(train_df.mean())
    corr_np = np.corrcoef(train_fill.values, rowvar=False)

    traits = train_df.columns.tolist()
    trait_index = {t: i for i, t in enumerate(traits)}

    print("# building predictor map")
    target_traits = df.columns[df.isna().sum() > 0]

    predictor_map, mean_cor_map = build_predictor_map(traits,trait_index,corr_np,top_k,cor_threshold)

    total = len(target_traits)
    print(f"# traits with missing: {total}")

    sub_train_map = {}
    for trait in target_traits:
        predictors = predictor_map[trait]
        if len(predictors) == 0:
            sub_train_map[trait] = train_df[[trait]]
        else:
            sub_train_map[trait] = train_df[predictors + [trait]]

    results = Parallel(n_jobs=threads)(
        delayed(impute_trait)(sub_train_map[trait],mean_cor_map[trait],trait,n_iter,i + 1,total,method)
        for i, trait in enumerate(target_traits)
    )

    df_out = df.copy()
    for trait, series in results:
        df_out[trait] = series

    df_out.to_csv(output, sep="\t")

    elapsed = time.time() - start_time
    print(
        f"# total time: {int(elapsed//3600):02d}:"
        f"{int(elapsed%3600//60):02d}:"
        f"{int(elapsed%60):02d}"
    )

    return df_out


def iterative_rf_single(sub_df, target_trait, n_iter=5,method="iterative_rf"):
    df = sub_df.copy()
    missing_mask = df[target_trait].isna()
    if missing_mask.sum() == 0:
        return df[target_trait]
    predictors = df.columns.drop(target_trait)
    df[predictors] = df[predictors].fillna(df[predictors].mean())
    df[predictors] = df[predictors].replace([np.inf, -np.inf], np.nan).fillna(0)
    df[target_trait] = df[target_trait].fillna(df[target_trait].mean())
    for _ in range(n_iter):
        train = df.loc[~missing_mask]
        test = df.loc[missing_mask]
        X_train, y_train = train[predictors], train[target_trait]
        X_test = test[predictors]
        if method == "iterative_rf":
            model = RandomForestRegressor(n_estimators=100, n_jobs=1, random_state=0)        
        elif method == "iterative_bayes":
            model = BayesianRidge()
        #model = RandomForestRegressor(n_estimators=100, n_jobs=1, random_state=0)
        model.fit(X_train, y_train)
        df.loc[missing_mask, target_trait] = model.predict(X_test)
    return df[target_trait]


def build_predictor_map(traits, trait_index, corr_np, top_k, cor_threshold):
    predictor_map = {}
    mean_cor_map = {}
    for trait in traits:
        i = trait_index[trait]
        cor = corr_np[i]
        idx = np.where(np.abs(cor) > cor_threshold)[0]
        idx = idx[idx != i]
        if len(idx) == 0:
            predictor_map[trait] = []
            mean_cor_map[trait] = 0
            continue
        cor_vals = np.abs(cor[idx])
        order = np.argsort(-cor_vals)
        idx = idx[order][:top_k]
        predictors = [traits[j] for j in idx]
        predictor_map[trait] = predictors
        mean_cor_map[trait] = cor_vals[order][:top_k].mean()
    return predictor_map, mean_cor_map


def impute_trait(sub_train_df, mean_cor, trait, n_iter, idx, total,method="iterative_rf"):
    missing_mask = sub_train_df[trait].isna()
    missing_n = missing_mask.sum()
    predictors = sub_train_df.columns.drop(trait)
    if len(predictors) == 0:
        print(f"# {idx}/{total} trait={trait} predictors=0 missing={missing_n} mean_cor=0")
        return trait, sub_train_df[trait]
    imputed = iterative_rf_single(sub_train_df, trait, n_iter,method)
    print(f"# {idx}/{total} trait={trait} predictors={len(predictors)} "
          f"missing={missing_n} mean_cor={mean_cor:.3f}")
    return trait, imputed

#----- -----  ----- ----- ----- ----- cluster ----- ----- ----- ----- ----- -----     
def compute_pca(data):
    # case 1: only one trait
    if data.shape[1] < 2:
        pc1 = data.iloc[:, 0].values.astype(float)

        # standardize PC1
        pc1 = (pc1 - pc1.mean()) / pc1.std(ddof=0)

        pc1_df = pd.DataFrame(
            pc1,
            index=data.index,
            columns=["PC1"]
        )

        return 1.0, pc1_df

    # case 2: multiple traits
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data)

    pca = PCA(n_components=1)
    pc1_scores = pca.fit_transform(data_scaled)[:, 0]

    explained_var = pca.explained_variance_ratio_[0]

    # standardize PC1
    pc1_scores = (pc1_scores - pc1_scores.mean()) / pc1_scores.std(ddof=0)

    pc1_df = pd.DataFrame(
        pc1_scores,
        index=data.index,
        columns=["PC1"]
    )

    return explained_var, pc1_df


def plot_group_bar(file_path, output_pdf):
    df = pd.read_csv(file_path, sep="\t", header=None, engine='python')
    
    df = df.iloc[:20, :3]
    df.columns = ["Group", "Count", "Percentage"]
    
    df["Percentage"] = df["Percentage"].str.rstrip('%').astype(float)
    df["Count"] = df["Count"].astype(int)
    
    df = df.sort_values("Count", ascending=False)
    
    plt.figure(figsize=(12, 6))
    
    sns.set_style("white")
    
    norm_perc = 0.2 + 0.8 * (df["Percentage"] - df["Percentage"].min()) / (df["Percentage"].max() - df["Percentage"].min())
    colors = sns.color_palette("Blues", n_colors=100)
    color_idx = (norm_perc * 99).astype(int)
    bar_colors = [colors[i] for i in color_idx]
    
    bars = plt.barh(df["Group"], df["Count"], color=bar_colors)
    
    for bar, perc in zip(bars, df["Percentage"]):
        plt.text(bar.get_width() + max(df["Count"])*0.01, 
                 bar.get_y() + bar.get_height()/2,
                 f"{perc:.2f}%", va='center', fontsize=10)
    
    plt.gca().invert_yaxis()    
    plt.xlabel("Count")
    plt.ylabel("Group")
    plt.title("Group Counts with PC1")
  
    plt.tight_layout()
    plt.savefig(output_pdf)
    plt.close()
    
    print(f"Bar plot saved to {output_pdf}")


def cor_fig(df,save_path):
    correlation_matrix = df.corr()
    plt.figure(figsize=(10,10))
    sns.clustermap(correlation_matrix,cmap="coolwarm", figsize=(10, 10))
    plt.savefig(f"{save_path}.png", dpi=300)
    correlation_matrix.to_csv(f"{save_path}.cor_matrix.tsv",sep="\t")



if __name__ == "__main__":
    main()
