#!/usr/bin/env python3 
import re
import sys 
import logging
import argparse
import numpy as np
import pandas as pd
import Common
import Evidence
import matplotlib.pyplot as plt
from scipy import stats
#import pysam
from cyvcf2 import VCF


def get_args(args):
    parser = argparse.ArgumentParser(description="Processing muti-traits GWAS with Emmax form vcf file",
                                    formatter_class=Common.CustomFormatter
        )
    parser.add_argument("-i", "--input",required=True, help="Path of the Gwas or Other site-value Result SNPs, such as *emmax.ps.tsv ",type=str)
    parser.add_argument("-g", "--gene",required=True, help="candidate gene name",type=str)
    parser.add_argument("-p", "--phe",required=True, help="Path of the input phenotype file",type=str)
    parser.add_argument("-v", "--vcf",required=True, help='vcf file in bgzip format and has bcftools-csi index',type=str)
    parser.add_argument("-a", "--annovar",required=True,help='Prefix of snp Annovar result file, *variant_function and *exonic_variant_function',type=str)
    parser.add_argument("-b","--bed",required=True,help="Table split bed file of gene or mRNA with at least 5 column: Chrom start end strand gene_name",type=str)
    parser.add_argument("-r","--region",help='Region(s) to run Haplotype, can be multiple,choices: All,CDS,NoSyn,intronic,Promoter ',type=str,choices=["All","CDS","NoSyn","intronic","Promoter"],nargs="+",default=["All"])
    parser.add_argument("-m","--method",help='Method to run Haplotype',type=str,choices=["Hom","Phase","Het"],default="Hom")
    parser.add_argument("-f","--flank",help='Flank size of Promoter(bp)',type=int,default=2000)
    parser.add_argument("-o", "--out",help="Prefix of outfiles",type=str,default="Hap")
    parser.add_argument("--chrom",help="The column of Chromosom name in the inputfile",type=str,default="Chrom") 
    parser.add_argument("--pos",help="The column of Maker pos in the inputfile",type=str,default="Pos")
    parser.add_argument("--value",help="The column of Value in the inputfile",type=str,default="Value")
    parser.add_argument("--name",help="The column of SNP id in the inputfile",type=str,default="ID")
    
    parsed_args = parser.parse_args(args)

    Common.check_path_exists(parsed_args.input)
    Common.check_path_exists(parsed_args.phe)
    Common.check_path_exists(parsed_args.vcf)
    Common.check_path_exists(parsed_args.bed)
    Common.check_path_exists(f"{parsed_args.annovar}.variant_function")
    Common.check_path_exists(f"{parsed_args.annovar}.exonic_variant_function")


    return parsed_args


def main(args=None):
    args=get_args(args)
    print (args)
    out= args.out


    print(f"Reading GWAS result from {args.input}")
    df=pd.read_table(args.input, sep="\t")
    name_col= {args.chrom:"CHROM",args.pos:"POS",args.value:"P",args.name:"ID"}
    snp_df = df.rename(columns=name_col)

    df_gene,df_prot = get_snps(snp_df,args.bed,args.gene,args.flank)

    (snp_bases,snp_pos2id) = Evidence.significance(df_gene)
    
    #print(snp_pos2id)
    print(f"Reading Annovar result")
    snp_genes = Evidence.read_annovar(f"{args.annovar}.variant_function",f"{args.annovar}.exonic_variant_function",snp_pos2id)
    df_snp_genes = pd.DataFrame(snp_genes,columns=["ID", "Gene", "Type", "Distance"]) 
    df_gene = df_gene.merge(df_snp_genes,on="ID",how="left")

    df_all = pd.concat([df_gene, df_prot],axis=0,ignore_index=True)
    df_all = df_all.sort_values("POS").reset_index(drop=True)
    print(df_all)

    df_candi=select_region(df_all,set(args.region))
    print("#Select SNP-info summary:")
    print(df_candi)

    sample_hap,sample_trait = haplotype_phenotype_analysis(df_candi,vcf_file=args.vcf, pheno_file=args.phe)

    hap_info_df, sample_hap_df = summarize_haplotypes(sample_hap, sample_trait,df_candi)

    df_all.to_csv(f"{out}.Site_info.txt",sep="\t",index=False)
    hap_info_df.to_csv(f"{out}.hap_info.txt",sep="\t",index=False)
    sample_hap_df.to_csv(f"{out}.sample_hap.txt",sep="\t",index=False)
    print(hap_info_df)
    print(sample_hap_df)



def select_region(df,region_set):

    if "All" in region_set:
        return df

    else:

        type_map = {
            "CDS": ["NoSyn", "Syn"],
            "NoSyn": ["NoSyn"],
            "Syn": ["Syn"],
            "Promoter": ["Promoter"],
            "Intron": ["intronic"]
        }
    
        selected_types = []
    
        for r in region_set:
            if r in type_map:
                selected_types.extend(type_map[r])
    
        df_candi = df[df["Type"].isin(selected_types)]
        return df_candi



def haplotype_phenotype_analysis(df_nosyn, vcf_file, pheno_file):
    """
    Construct gene-level haplotypes using homozygous nonsynonymous SNPs.
    Optimized: fetch each SNP by genomic position (region) instead of traversing entire VCF.

    Parameters
    ----------
    df_nosyn : pd.DataFrame
        Columns: CHROM, POS, ID, P, Gene, Type, Distance
    vcf_file : str
        Path to bgzipped & tabix-indexed VCF
    pheno_file : str
        Phenotype file: FID IID Trait

    Returns
    -------
    sample_hap : dict
        sample -> haplotype string
    sample_trait : dict
        sample -> phenotype
    """

    # -------------------------------
    # Step 0: load phenotype
    # -------------------------------
    pheno_df = pd.read_csv(pheno_file, sep=None, engine='python', header=None, names=["FID","IID","Trait"])
    sample_trait = dict(zip(pheno_df["IID"], pheno_df["Trait"]))

    # -------------------------------
    # Step 1: sort SNPs by position
    # -------------------------------
    df_nosyn = df_nosyn.sort_values("POS")
    vcf = VCF(vcf_file)
    samples = vcf.samples
    hap_temp = {s: [] for s in samples}

    # -------------------------------
    # Step 2: iterate each SNP by region fetch
    # -------------------------------
    for idx, row in df_nosyn.iterrows():
        chrom = row["CHROM"]
        pos = row["POS"]
        snp_id = row["ID"]

        # fetch 是 0-based start, 1-based end
        for variant in vcf('%s:%s-%s' %(chrom, pos, pos)):
            if variant.ID != snp_id:
                continue  # 确保 ID 对应
            for i, gt in enumerate(variant.genotypes):   #v.gt_types
                sample = samples[i]
                #print(gt)
                a1, a2, _ = gt
                if a1 == a2:
                    if a1 == 0:
                        hap_temp[sample].append(variant.REF)
                    else:
                        hap_temp[sample].append(variant.ALT[0])
                else:
                    hap_temp[sample].append("N")  # 杂合标记 N

    # -------------------------------
    # Step 3: 构建 haplotype
    # -------------------------------
    sample_hap = {}
    for sample, alleles in hap_temp.items():
        if "N" in alleles:
            continue  # 含杂合丢掉
        if sample in sample_trait:
            sample_hap[sample] = "|".join(alleles)

    print(sample_trait)
    return sample_hap,sample_trait


def summarize_haplotypes(sample_hap, sample_trait, df_nosyn):
    """
    Generate haplotype-level and sample-level summaries.

    Parameters
    ----------
    sample_hap : dict
        sample -> haplotype string (alleles joined by "|")
    sample_trait : dict
        sample -> phenotype
    df_nosyn : pd.DataFrame
        SNP annotation table with ID column

    Returns
    -------
    hap_info_df : pd.DataFrame
        Columns: HapID, Count, HapType, Site-info, Phe
    sample_hap_df : pd.DataFrame
        Columns: Sample, HapID, Phe
    """

    # -------------------------------
    # Step 1: 统计 haplotype 样本数
    # -------------------------------
    hap_to_samples = {}
    for sample, hap_str in sample_hap.items():
        if hap_str not in hap_to_samples:
            hap_to_samples[hap_str] = []
        hap_to_samples[hap_str].append(sample)

    # -------------------------------
    # Step 2: 按样本数排序生成 HapID
    # -------------------------------
    hap_sorted = sorted(hap_to_samples.items(), key=lambda x: len(x[1]), reverse=True)
    hap_info_list = []
    hapid_map = {}  # hap_str -> HapID

    for i, (hap_str, samples) in enumerate(hap_sorted):
        hap_id = f"hap{i+1}"
        hapid_map[hap_str] = hap_id
        phe_values = [sample_trait[s] for s in samples if s in sample_trait]
        mean_phe = sum(phe_values)/len(phe_values) if phe_values else None
        # Site-info
        site_info = ",".join(df_nosyn["ID"].tolist())  # 所有 SNP
        hap_info_list.append({
            "HapID": hap_id,
            "Count": len(samples),
            "HapType": hap_str,
            "Site-info": site_info,
            "Phe": mean_phe
        })

    hap_info_df = pd.DataFrame(hap_info_list)

    # -------------------------------
    # Step 3: 构建 sample-hap 表
    # -------------------------------
    sample_hap_list = []
    for sample, hap_str in sample_hap.items():
        hap_id = hapid_map[hap_str]
        phe = sample_trait.get(sample, None)
        sample_hap_list.append({
            "Sample": sample,
            "HapID": hap_id,
            "Phe": phe
        })

    sample_hap_df = pd.DataFrame(sample_hap_list)

    return hap_info_df, sample_hap_df



def get_snps(df,bedfile,gene,flank):
    gene_location = Common.read_file(bedfile,mode="dict",keys=[4],vals=[0,1,2,3])

    if gene not in gene_location:
        print(f"# Error: {gene} not in bedfile: {bedfile}, please check it.")

    c,s,e,strand =  gene_location[gene]

    df_gene = df[(df["CHROM"] == c ) & (df["POS"] >= int(s)) &  (df["POS"] <= int(e)) ].copy()

    if  strand  == "+":
        df_prot = df[(df["CHROM"] == c ) & (df["POS"] > int(s) - flank ) &  (df["POS"] < int(s)) ].copy()
        df_prot["Gene"] = gene
        df_prot["Type"] = "Promoter"
        df_prot["Distance"] = int(s) - df_prot["POS"] + 1
    elif strand  == "-":
        df_prot = df[(df["CHROM"] == c ) & (df["POS"] > int(e)) &  (df["POS"] < int(e) + flank) ].copy()
        df_prot["Gene"] = gene
        df_prot["Type"] = "Promoter"
        df_prot["Distance"] =  df_prot["POS"] - int(e) + 1
    else:
        print(f"None strand info of gene in bedfile {bedfile}")

    return df_gene,df_prot


if __name__ == '__main__':
    #logger = configure_logging()
    main()
