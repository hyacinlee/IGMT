#!/usr/bin/env python3 
import re
import sys 
import logging
import argparse
import numpy as np
import pandas as pd
import Common
import statsmodels.api as sm
import matplotlib.pyplot as plt
#from qmplot import manhattanplot
from scipy import stats
from cyvcf2 import VCF

#logger=""

def get_args(args):
    parser = argparse.ArgumentParser(description="Processing muti-traits GWAS with Emmax form vcf file",
                                    formatter_class=Common.CustomFormatter
        )
    parser.add_argument("-i", "--input",required=True, help="Path of the Gwas or Other site-value Result of Candidata SNPs, such as *significant.tsv ",type=str)
    parser.add_argument("-v", "--vcf",required=True, help='vcf file in bgzip format and has bcftools-csi index',type=str)
    parser.add_argument("-a", "--annovar",required=True,help='Prefix of snp Annovar result file, *variant_function and *exonic_variant_function',type=str)
    parser.add_argument("-b","--bed",required=True,help="Table split bed file of gene or mRNA with at least 5 column: Chrom start end strand gene_name",type=str)
    parser.add_argument("-p", "--phe", help="Path of the input phenotype file",type=str)
    parser.add_argument("-f", "--function",help='Function Annotation of Genes. Support multi-column information, but the first column must be the gene name and include the table header',type=str)
    parser.add_argument("-e", "--eQTL",help='Known eQTLs results file in SNPs -- Gens formats',type=str)
    parser.add_argument("-o", "--out",help="Prefix of outfiles",type=str,default="Evi")
    parser.add_argument("-w", "--weight",help='Other evidence with weight',type=str)
    parser.add_argument("--chrom",help="The column of Chromosom name in the inputfile",type=str,default="Chrom") 
    parser.add_argument("--pos",help="The column of Maker pos in the inputfile",type=str,default="Pos")
    parser.add_argument("--value",help="The column of Value in the inputfile",type=str,default="Value")
    parser.add_argument("--name",help="The column of SNP id in the inputfile",type=str,default="ID")
    parser.add_argument("--min_sites",help="The minimum number of Candidata SNPs in the input file",type=int,default=5)
    parser.add_argument("--max_sites", help="The maximum number of candidate SNPs used for Evidence analysis",type=int,default=200)
    parser.add_argument("--eQTL_score",help="The weights score of eQTL SNPs",type=int,default=3)
    parser.add_argument("--regions",help="The range size between the leader-SNPs and candidate genes, in Kb",nargs="+",default=[100])
    parser.add_argument("--regions_score",help="The weights scores of far region SNPs,must be consistent with the number of the --regions",nargs="+",default=[1])
    parser.add_argument("--snp_scores",help="The weights of the following various type SNPs: [NoSys, Sys, intronic , splicing]",nargs="+",default=[10,5,3,2])
    parser.add_argument("--base_scores",help="The weights of the Basic GWAS result (range from 0~1)",type=float,default=0.8)
    

    parsed_args = parser.parse_args(args)

    Common.check_path_exists(parsed_args.input)
    #Common.check_path_exists(parsed_args.phe) # Adding the support for phe lacking. 20260302
    Common.check_path_exists(parsed_args.vcf)
    Common.check_path_exists(parsed_args.bed)
    Common.check_path_exists(f"{parsed_args.annovar}.variant_function")
    Common.check_path_exists(f"{parsed_args.annovar}.exonic_variant_function")

    if not len(parsed_args.regions) == len(parsed_args.regions_score):
        print("# Error: --regions_score must hava the same length with --regions")
        exit(1)


    if parsed_args.function:
        Common.check_path_exists(parsed_args.function)    

    return parsed_args


def main(args=None):
    args=get_args(args)
    print (args)
    out= args.out

    #logger = configure_logging()
    #logger.info('Starting main Program')
    #logger.info(args)
 
    #logger.info(f"Reading GWAS result from {args.input}")
    #print(f"Reading GWAS result from {args.input}")
    #df=pd.read_table(args.input, sep="\t")
    #if len(df) < args.min_sites:
    #    print(f"# Warning : Significance locus in file {args.input} was less than {args.min_sites}, so this program is suspended here")
    #    exit(0)

    #name_col= {args.chrom:"CHROM",args.pos:"POS",args.value:"P",args.name:"ID"}
    #snp_df = df.rename(columns=name_col)
    
    snp_df = judge_df(args)

    gene_header,base_score_dict = check_header(args.regions,args.regions_score,args.eQTL,args.eQTL_score,args.snp_scores)

    (snp_bases,snp_pos2id) = significance(snp_df)
    
    #logger.info(f"Reading Candidata SNP Genotype info from {args.vcf}")
    print(f"Reading Candidata SNP Genotype info from {args.vcf}")
    (snp_GTs,snp_bases) = read_GTs(args.vcf,snp_bases,out)

    print(f"Reading Annovar result")
    snp_genes = read_annovar(f"{args.annovar}.variant_function",f"{args.annovar}.exonic_variant_function",snp_pos2id)

    print(f"Reading Gene feature bed to identify flank genes")
    snp_genes = annot_genes(args.bed,args.regions,snp_pos2id,snp_genes)

    if args.eQTL:
        snp_genes=annot_genes_eqtl(args.eQTL,snp_bases,snp_genes,args.bed)

    #print(snp_genes)
    #logger.info("Reading Phenotype data from %s" %(args.phe)) 
    if args.phe:
        Phe = Common.read_file(args.phe,mode="dict",keys=[0],vals=[2])
    else:
        Phe = None # Adding the support for phe lacking. 20260302
    #Phe = readTableDict(args.phe,0,2)

    (gene_snp_count,gene_snp_count_cor)=out_sites_genes(out,snp_bases,snp_genes,snp_GTs,Phe,snp_df)

    df_base = out_basic_genes(out,gene_snp_count,gene_snp_count_cor,gene_header,base_score_dict,args.function)

    if args.weight:
        out_weight_genes(out, df_base, args.weight, args.base_scores, Phe, args.bed, args.function)


def judge_df(args):
    """
    Read GWAS results, check the number of significant SNPs,
    and retain at most max_site SNPs ranked by P-value.

    Returns
    -------
    snp_df : pandas.DataFrame
        Filtered candidate SNP dataframe.
    """
    input_file = args.input
    out        = args.out
    min_site   = args.min_sites
    max_site   = args.max_sites

    print(f"Reading GWAS result from {input_file}")
    df = pd.read_table(input_file, sep="\t")

    # Check minimum number of candidate SNPs
    if len(df) < min_site:
        print(f"# Warning : Significant loci in file {input_file} were fewer than {min_site}, so this program is suspended here")
        exit(0)

    # Standardize column names
    name_col = {
        args.chrom: "CHROM",
        args.pos: "POS",
        args.value: "P",
        args.name: "ID"
    }
    snp_df = df.rename(columns=name_col)

    # Make sure P-value is numeric
    snp_df["P"] = pd.to_numeric(snp_df["P"], errors="coerce")
    snp_df = snp_df.dropna(subset=["P"])

    # Sort all candidate SNPs by P-value
    snp_df = snp_df.sort_values(
        by="P",
        ascending=True
    ).reset_index(drop=True)

    # Retain only the top N SNPs if necessary
    if max_site > 0 and len(snp_df) > max_site:
        print(
            f"# Info : {len(snp_df)} candidate SNPs were detected. "
            f"Only the top {max_site} SNPs ranked by P-value "
            f"will be used for Evidence analysis."
        )
        snp_df = snp_df.head(max_site).copy()

    # Output SNPs actually used for Evidence analysis
    valid_file = f"{out}.sign.vaildSite.xls"
    snp_df.to_csv(
        valid_file,
        sep="\t",
        index=False
    )

    print(f"# Info : Valid SNPs were written to {valid_file}")
    return snp_df



def out_weight_genes(out, df, weight, base_scores, Phe, bed, function):
    
    evidens=Common.read_file(weight,mode="list",vals=[0,1,2,3])

    print(evidens)
    #print(df)
    second_score_dict={}
    for name,model,value,file in evidens:
        print(f"Evidence gene from {model} file {file}")
        if model == "express":
            df=evidence_express(file,df,name,Phe)     # Phe vs FPKM
            #second_score_dict.append([name,value])      
        elif model == "gene":   # Known Candidate Gene
            df=evidence_genes(file,df,name)
            #second_score_dict.append([name,value])          
        elif model == "loci":   # Known Candidate locis 
            df=evidence_loci(file,df,name,bed)
            #second_score_dict[name]=value
        else:
            print(f"# Unknow evidence data type {model},Ignore and skip file {file}") 
        second_score_dict[f"{name}_v"]=value

    
    df = cal_df_score(df,second_score_dict,"Evi_score")
    add_dict=dict(zip(["Base_score","Evi_score"],[base_scores,1-base_scores]))
       
    df = cal_df_score(df,add_dict,"Total_score")

    df = re_order_df(df, ['Total_score', 'Base_score', 'Evi_score','Total_SNP'])
    df.sort_values(by='Total_score', ascending=False, inplace=True)

    for key in second_score_dict.keys():
        df = df.drop(key, axis=1)

    if function:
        df=add_function(df,function)
    print(df)
    df.to_csv(f"{out}.sign.WeightGenes.xls",sep="\t")    


def out_sites_genes(out,snp_bases,snp_genes,snp_GTs,Phe,snp_df):
    gene_snp_count={}
    gene_snp_count_cor={}
    #   print(snp_df)
    #logger.info("Outputing significance SNP function, GWAS-Pvalue")
    with open("%s.sign.BasicSite.xls" % (out),"w") as of1:
        if "R2" in snp_df.columns:
            of1.write("ID\tChr\tPos\tRef\tAlt\tGene\tType\tDistance\tR2\tP-GWAS\tP-CorPhe\n")
        else:
            of1.write("ID\tChr\tPos\tRef\tAlt\tGene\tType\tDistance\tP-GWAS\tP-CorPhe\n")

        for rs in sorted(snp_genes,key=lambda x: snp_bases[x[0]][2]):
            #    #logger.warning("Check SNP %s %s" %(rs[0],rs[1]))
            sID  = rs[0]
            base = snp_bases[sID]
            #print(snp_GTs[sID],Phe)
            if Phe:
                (pp1,p1,r)=stTest(snp_GTs[sID],Phe)     # student t-test SNP vs Phe
            else:
                (pp1,p1,r)=(0,0,1)    # Adding the support for phe lacking. 20260302

            if "R2" in snp_df.columns:
                r2 = snp_df.loc[snp_df['ID'] == sID, 'R2'].values[0].round(4)
                of1.write(f"{rs[0]}\t{base[0]}\t{base[1]}\t{base[3]}\t{base[4]}\t{rs[1]}\t{rs[2]}\t{rs[3]}\t{r2}\t{base[2]}\t{p1}\n")
            else:
                of1.write(f"{rs[0]}\t{base[0]}\t{base[1]}\t{base[3]}\t{base[4]}\t{rs[1]}\t{rs[2]}\t{rs[3]}\t{base[2]}\t{p1}\n")

            gene_snp_count=Common.accumulateDict(gene_snp_count,1,rs[1],rs[2])
            gene_snp_count=Common.accumulateDict(gene_snp_count,1,rs[1],"Total_SNP")
            gene_snp_count_cor=Common.accumulateDict(gene_snp_count_cor,0,rs[1],rs[2])
            gene_snp_count_cor=Common.accumulateDict(gene_snp_count_cor,0,rs[1],"Total_SNP")
            if p1 < 0.01:
                gene_snp_count_cor=Common.accumulateDict(gene_snp_count_cor,1,rs[1],rs[2])
                gene_snp_count_cor=Common.accumulateDict(gene_snp_count_cor,1,rs[1],"Total_SNP")

    return gene_snp_count,gene_snp_count_cor
    

def out_basic_genes(out,gene_snp_count,gene_snp_count_cor,gene_header,base_score_dict,function):

    gene_snp_count = Common.fill_double_dict_value(gene_snp_count,gene_header,0)

    gene_snp_count_cor = Common.fill_double_dict_value(gene_snp_count_cor,gene_header,0)

    out_dict={}
    for gene in gene_snp_count.keys():
        outs = [ gene_snp_count[gene][x] for x in gene_header]
        out_dict[gene] = dict(zip(gene_header,outs))
    
    df_out = pd.DataFrame(out_dict).T
    df_out.index.name="Gene"

    leaders=leader_snp(f"{out}.sign.BasicSite.xls")

    df_cal = cal_df_score(pd.DataFrame(gene_snp_count_cor).T,base_score_dict,"Base_score")
    df_cal.index.name="Gene"
    #print(df_cal)

    df_out["Base_score"]=df_cal["Base_score"]

    df_out=df_out.join(leaders,how='left')
    #print(df_out)

    df_out.sort_values(by='Base_score', ascending=False, inplace=True)
    df_nofunction = df_out.copy()

    df_out = re_order_df(df_out, ['Base_score','Total_SNP'])

    if function:
        df_out=add_function(df_out,function) 

    df_out.to_csv(f"{out}.sign.BasicGene.xls",sep="\t")

    return df_nofunction   



def leader_snp(dfile):
    df = pd.read_csv(dfile, sep="\t")

    #print(df)
    def type_priority(t):
        t = t.lower()
        if t == "syn":
            return 0
        elif t == "nosyn":
            return 1
        elif t == "intronic":
            return 2
        elif t == "splicing":
            return 3
        elif t == "eqtl":
            return 4
        else:
            match = re.search(r"(\d+)kb", t)
            return 5 + (int(match.group(1)) if match else 10000)

    df["TypePriority"] = df["Type"].apply(type_priority)

    results = []

    for gene, subdf in df.groupby("Gene"):
        # Highest_SNP: P-GWAS最小值对应的记录
        highest_row = subdf.loc[subdf["P-GWAS"].idxmin()]
    
        # Nearest_SNP: 多重排序（Type + P-GWAS）
        sorted_subdf = subdf.sort_values(by=["TypePriority", "P-GWAS"])
        nearest_row = sorted_subdf.iloc[0]
    
        results.append({
            "Gene": gene,
            "Highest_ID": highest_row["ID"],
            "Highest_Type": highest_row["Type"],
            "Highest_P": highest_row["P-GWAS"],
            "Highest_Dis": highest_row["Distance"],
            "Nearest_ID": nearest_row["ID"],
            "Nearest_Type": nearest_row["Type"],
            "Nearest_P": nearest_row["P-GWAS"],
            "Nearest_Dis": nearest_row["Distance"]
        })
    
    # 转为DataFrame查看或保存
    leaders = pd.DataFrame(results)
    leaders = leaders.set_index("Gene")
    leaders.index.name="Gene"
    return(leaders)


def re_order_df(df, cols_to_front):
    new_cols = cols_to_front + [col for col in df.columns if col not in cols_to_front]
    df = df[new_cols]
    return df  


def annot_genes_eqtl(eQTL_file,snp_bases,snp_genes,bedfile):

    print(f"# eQTL evidence founded, use eQTL file: {eQTL_file} to identifed candidate genes")
    Common.check_path_exists(eQTL_file)
    eqlt_snps  = Common.read_file_accumulateDict(eQTL_file,vals=[1],key1=[0])

    gene_location = Common.read_file(bedfile,mode="dict",keys=[4],vals=[0,1,2,3])

    for ss in eqlt_snps.keys():
        if ss in snp_bases:
            for g in eqlt_snps[ss]:
                if g not in gene_location:
                    dis = "None"
                else:
                    s=int(gene_location[g][1])
                    e=int(gene_location[g][2])
                #print(snp_bases[ss])
                    cid,pos = snp_bases[ss][0],snp_bases[ss][1]
                    if not cid  == gene_location[g][0]:
                        dis = "Trans"
                    else:
                        if int(pos)  > e:
                            dis = int(pos)-e+1  
                        elif int(pos) < s:
                            dis = s - int(pos)+1
                        else:
                            dis=0
                snp_genes.append([ss,g,"eQTL",dis])

    return snp_genes


def check_header(regions,regions_score,eQTL,eQTL_score,snp_scores):
    gene_header=["NoSyn","Syn","intronic","splicing"]
    for i in range(len(regions)):
        gene_header.append(f"{regions[i]}Kb")
        snp_scores.append(int(regions_score[i]))
    if eQTL:
        gene_header.append("eQTL")
        snp_scores.append(int(eQTL_score))

    base_score_dict = dict(zip(gene_header,snp_scores))
    gene_header.insert(0,"Total_SNP")
    return gene_header,base_score_dict



def add_function(df,function_file):
    functions,header=Common.read_file(function_file,mode="dict",vals=[],keys=[0],header=True,sep="\t")
    for i in range(1,len(header)):
        info_dict={}
        for g in functions.keys():
            info_dict[g] = functions[g][i]
        df=updata_df(df,header[i],info_dict,missing="---")

    return df


def cal_df_score(df,weights,name):
    #print(df)
    #print(weights)
    weight_dict = {col: float(val) for col, val in weights.items()}
    df[f'{name}_tmp'] = sum(df[col] * weight for col, weight in weight_dict.items())
    df[f'{name}'] = df[f'{name}_tmp'] / df[f'{name}_tmp'].sum() *100
    df[f'{name}'] = df[f'{name}'].round(4)
    df = df.drop(f'{name}_tmp', axis=1)
    return df


def evidence_loci(file,df_handle,name,bedfile):
    if not bedfile:
        print("# Error: Please specify the bed file containing mRNA or genes using the --bed !")
        exit(1)

    gene_location = Common.read_file(bedfile,mode="dict",keys=[4],vals=[0,1,2,3])
    #print(gene_loci)
    #locis = Common.read_file(file,mode="dict",vals=[0,1,2,3],keys=[0])   # chrom start end id
    locis = Common.read_file_accumulateDict(file,vals=[0,1,2,3],key1=[0])

    gene_value={}
    gene_info={}
    candidate_genes=df_handle.index.tolist()
    for gene in candidate_genes:
        #print(gene)
        chrid=gene_location[gene][0]
        overlapping_regions=Common.search_cloesd_region(gene_location[gene],locis[chrid])
        if len(overlapping_regions) >= 1:
            gene_value[gene] = len(overlapping_regions)
            gene_info[gene] = ",".join([x[3] for x in overlapping_regions])
        else:
            gene_value[gene] = 0
            gene_info[gene] = "None"
        #print(gene,overlapping_regions)
    df_handle=updata_df(df_handle,f"{name}_v",gene_value)
    df_handle=updata_df(df_handle,name,gene_info)
    return df_handle


def evidence_genes(file,df_handle,name):
    info={}
    value={}
    try:
        genes_evi = Common.read_file(file,mode="dict",keys=[0],vals=[0,1,2],sep="\t")
    except Exception:
        print(f"###")
        print(f"# Reading Error! {file} must have 3 columns: gene_id value info")
        print(f"###")
        exit(1)

    for gene in df_handle.index.tolist():
        if gene in genes_evi:
            value[gene] = float(genes_evi[gene][1])
            info[gene]  = genes_evi[gene][2]
        else:
            value[gene] = 0
            info[gene]  = "None"

    df_handle=updata_df(df_handle,f"{name}_v",value)
    df_handle=updata_df(df_handle,name,info)
    return df_handle


def evidence_express(file,df_handle,name,Phe):
    #logger.info("Reading the gene expression matrix from %s" %(args.express))
    gene_asso={}
    gene_value={}
    G,Gsample,Genes  = readData(file,"v")
    candidate_genes = df_handle.index.tolist()
    for gene in candidate_genes:
        gene_value[gene] = 0
        if gene not in G:
            print("%s has no Transcript Data, skipping" % (gene))
            gene_asso[gene] = "NA"
        else:
            lineregress=expressPheAsso(G[gene],Phe)
            #print(gene,lineregress.pvalues,lineregress.pvalues[1])
            gene_value[gene] = 0
            gene_asso[gene]  ="None"
            if not np.isnan(lineregress.pvalues[1]) and lineregress.pvalues[1] < 0.01:
                gene_value[gene] = 1
                gene_asso[gene]  = lineregress.pvalues[1]
            #print(gene,lineregress.pvalues,lineregress.pvalues[1],p)
            #gene_asso[gene] = p

    df_handle=updata_df(df_handle,name,gene_asso)
    df_handle=updata_df(df_handle,f"{name}_v",gene_value)

    return df_handle


def updata_df(df,index_name,data_dict,missing="None"):
    df[index_name] = df.index.map(data_dict)
    df[index_name].fillna(missing, inplace=True)
    #new_order = [c for c in df.columns if c != 'B'] + ['B']
    return df



def expressPheAsso(expressDict,pheDict,gname=None,tname=None):
    x=[]
    y=[]
    for sam in pheDict.keys():
        if pheDict[sam] == "NA" or sam not in expressDict:
            continue
        else:
            x.append(float(expressDict[sam]))
            y.append(float(pheDict[sam]))
    x = sm.add_constant(x)
    model = sm.OLS(y,x).fit()

    if not gname == None:
        plt.style.use('seaborn')
        fig,ax = plt.subplots()
        ax.scatter(x, y)
        ax.set_xlabel(gname, fontsize=12)
        ax.set_ylabel(tname, fontsize=12)
        ax.tick_params(axis='both', which='major', labelsize=10)
        qqline(ax, "r", x, y)
        plt.savefig('%s_vs_%s.png' %(tname,gname))

    return model


def stTest(gt,Vals):
    gts={}
    gtl=[]
    for sample in gt.keys():
        if sample not in Vals:
            continue 
        if Vals[sample] == "NA":
            continue
        gtl=accumulateList(gtl,gt[sample])
        gts=accumulateDict(gts,str(Vals[sample]),gt[sample])

    gtl=sorted(gtl)
    p1,p2,p3 = (1,1,1)
    if len(gtl) == 3:
        t1,p1 = stats.ttest_ind(np.array(listType(gts[gtl[0]],"float")),np.array(listType(gts[gtl[1]],"float")))
        t2,p2 = stats.ttest_ind(np.array(listType(gts[gtl[1]],"float")),np.array(listType(gts[gtl[2]],"float")))
        t3,p3 = stats.ttest_ind(np.array(listType(gts[gtl[0]],"float")),np.array(listType(gts[gtl[2]],"float")))
    elif len(gtl) == 2:
        t1,p1 = stats.ttest_ind(np.array(listType(gts[gtl[0]],"float")),np.array(listType(gts[gtl[1]],"float")))
    
    p=min(p1,p2,p3)
    if p1 <0.01 or  p2 <0.01 or p3 <0.01:
        return 1,p,gts
    else:
        return 0,p,gts

# Reading Function

def significance(df):

    snp_bases={}
    snp_pos2id={}

    signDf=df
 
    for index,row in signDf.iterrows():
        c=str(row['CHROM'])
        p=int(row['POS'])
        i=str(row['ID'])
        v=float(row['P'])
        snp_bases[i]=[c,p,v]
        tags=f"{c}#{p}"
        snp_pos2id[tags] = i

    return (snp_bases,snp_pos2id)



def read_GTs(gzvcf,snp_bases,out=None):
    snp_GTs={}
    snpsINFO={}
    vcf_reader = VCF(gzvcf)
    samples=vcf_reader.samples
    with open("%s.sign.BasicGT.lst" %(out),"w") as of:
        of.write("ID\tChr\tPos\tRef\tAlt\t%s\n" %("\t".join(samples)))
        for k in snp_bases.keys():
            l=snp_bases[k]
            for v in vcf_reader('%s:%s-%s' %(l[0],l[1],l[1])):
                snp_GTs[k] =dict(zip(samples,v.gt_types))
                snp_bases[k].append(v.REF)
                snp_bases[k].append(v.ALT[0])
                types={0:f"{v.REF}{v.REF}",1:f"{v.REF}{v.ALT[0]}",2:f"NN",3:f"{v.ALT[0]}{v.ALT[0]}"}
                gts = [ types[x] for x in v.gt_types.tolist()]
                of.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(k,l[0],l[1],v.REF,v.ALT[0],"\t".join(gts)))  #gt_types(0 ref hom,1 het, 2 unknow ,3 alt hom) 

    return snp_GTs,snp_bases


def read_annovar(vf,evf,snp_pos2id):
    #print(snp_pos2id)
    snp_genes=[]
    snp_exon=[]
    for l in open(vf,"r"):
        ls =l.strip().split()

        if ls[0] =="intergenic" or "stream" in ls[0] or "UTR" in ls[0]:
            continue

        tags = "%s#%s" %(ls[2],ls[3])
        
        if ls[6] == "-": #  del, shift 1bp
            tags = f"{ls[2]}#{int(ls[3])-1}" 
        #if len(ls[4]) < len(ls[3]):  #  del, shift 1bp
        #    tags = f"{ls[2]}#{int(ls[3])-1}"
        # gene,loci    
        if tags not in snp_pos2id:
            #print(tags)
            continue

        #print(tags,l)
        pgs=re.split(r'\(|,',ls[1]) 

        if ls[0] == "exonic":
            for gene in pgs:
                snp_exon.append(tags)
        elif ls[0] == "splicing":
            snp_genes.append([snp_pos2id[tags],pgs[0],ls[0],0])
        else:                       
            for gene in pgs:
                snp_genes.append([snp_pos2id[tags],gene,ls[0],0])


    for le in open(evf,"r"):
        les=le.strip().split("\t")
        tags="%s#%s"%(les[3],les[4])
        #print (tag)
        if tags in snp_exon:
            #print(tags,les)
            gene=les[2].split(":")[0]
            if gene == "UNKNOW":
                print(f"# UNKNOW gene information in {evf}, please check it.")
                exit(1)
            #print (le)
            #if "nonsynonymous"  or "nonframeshift" in les[1]:  # indels trait as Syn
            if "nonsynonymous" in les[1]:
                snp_genes.append([snp_pos2id[tags],gene,"NoSyn",0])
            else:
                snp_genes.append([snp_pos2id[tags],gene,"Syn",0])

    return snp_genes


def annot_genes(bedfile,regions,snp_pos2id,snp_genes):

    if not bedfile:
        print("# Error: Please specify the bed file containing mRNA or genes using the --bed !")
        exit(1)

    gene_location = Common.read_file(bedfile,mode="dict",keys=[4],vals=[0,1,2,3])
    #print(gene_location)

    for gene in gene_location:
        s=int(gene_location[gene][1])
        e=int(gene_location[gene][2])
        for pos_inf in snp_pos2id:
            cid,pos = pos_inf.split("#")
            if gene_location[gene][0] == cid:
                for region_size in regions:
                    if int(pos) - e < int(region_size)*1000 and  int(pos) > e:    # s e pos  
                        snp_genes.append([snp_pos2id[pos_inf],gene,f"{region_size}Kb",int(pos)-e+1])
                    elif s - int(pos) < int(region_size)*1000 and int(pos) <s:    # pos s e 
                        snp_genes.append([snp_pos2id[pos_inf],gene,f"{region_size}Kb",s-int(pos)+1])
                    else:
                        continue

    return snp_genes



def readData(infile,order):
    P={}
    Hlist=[]
    Vlist=[]
    for l in open(infile,"r"):
        if "t_name" in l:
            continue

        ls = l.strip().split()
        name=ls.pop(0)

        if len(Hlist) ==0:
            Hlist=ls
        else:
            Vlist.append(name)
            for i in range(len(Hlist)):
                k1 = Hlist[i]
                k2 = name
                v = ls[i]
                if order == "h":  ##  H-line  as key
                    P=accumulateDict(P,v,k1,k2) 
                if order == "v":  ##  V-line  as key
                    P=accumulateDict(P,v,k2,k1)
    return P,Hlist,Vlist

def readTableDict(infile,ck=0,cv=1):
    info={}
    with open(infile,"r") as inf:
        for line in inf:
            ss = line.strip().split("\t")
            if cv != -1:
                info[ss[ck]] = ss[cv]
            else:
                info[ss[ck]] = line
    return info

# Gernal Function 
def listType(mylist,types):
    nl=[]
    for s in mylist:
        if types == "str":
            nl.append(str(s))
        elif types == "float":
            nl.append(float(s))
    return nl

def accumulateList(myList,val):
    if val not in myList:
        myList.append(val)
    return myList



def accumulateDict(myDict,val,keyA,keyB="no"):
    '''
        accumulate a 1~2d Dictionary
    '''
    if keyB == "no":
        if keyA in myDict:
            if type(val) == int :
                myDict[keyA] += val
            else:
                myDict[keyA].append(val)         

        else:
            if type(val) == int :
                myDict[keyA] = val
            else:
                myDict[keyA] = [val]

    else:
        if keyA in myDict:
            if keyB in myDict[keyA]:
                myDict[keyA][keyB] += val
            else:
                myDict[keyA].update({keyB: val})   
        else:
            myDict.update({keyA:{keyB: val}})

    return myDict



def configure_logging():
    #logging.basicConfig(filename="EGCV.log",level=logging.INFO,format='%(asctime)s %(levelname)s %(name)s %(message)s')
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(name)s %(message)s',
        datefmt="%Y-%m-%d %H:%M:%S",)
    return logging.getLogger("EVCG")


if __name__ == '__main__':
    #logger = configure_logging()
    main()
