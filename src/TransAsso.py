#!/usr/bin/env python3 
import sys 
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
#from statsmodels.graphics.gofplots import qqline
#from scipy.stats import pearsonr
import Common
import argparse

def get_args(args):
    parser = argparse.ArgumentParser(description="Association analysis between gene expression matrix and phenotype in population",
                                    formatter_class=Common.CustomFormatter
        )
    parser.add_argument("-e", "--express",required=True, type=str,help="Path of the input gene expression matrix file")
    parser.add_argument("-p", "--phe",required=True, type=str,help="Path of the input phenotype file, *phe")
    parser.add_argument("-c", "--cov",type=str,help="Covariance file")
    parser.add_argument("-t", "--trait",nargs="+",help="Trait Names to run TransAsso. By default, all traits in the Phe file")
    parser.add_argument("-l", "--loci",help="Path of the gene loci file in bed format with 5 Colums")
    parser.add_argument("-o", "--out",type=str,help="Out frefix",default="TransAsso")
    parser.add_argument("--threads",help="The number of parallel threads to perform GWAS",type=int,default=12)
    parser.add_argument("--figGenes",type=str, help="a list file, Out Fpkm - Phe correlation plot with figGenes list file")
    parsed_args = parser.parse_args(args)
    print (parsed_args)
    return parsed_args



def main(args=None):
    args=get_args(args)

    out=args.out    
    ### reads Phe 
    Phe = pd.read_table(args.phe,sep="\t",index_col=0)

    ### reads FPKM data
    Exp = pd.read_table(args.express,sep="\t",index_col=0)

    covar_df = read_Cov(args.cov)
    #print(covar_df.dtypes)
    #print(covar_df.head())


    Genes = Exp.index.tolist()
    print(f"## Total Genes: {len(Genes)}\n") 

    cut = 1.0/(len(Genes))
    with open(f"./Associate.cutline.txt","w") as cf:
        cf.write(str(cut))

    Traits = args.trait if args.trait else Phe.columns.tolist()

    tasks=[]
    for t in Traits:
        if t not in Phe.columns.tolist():
            print("# Warning: f{t} not match {args.phe} columns ")
            exit(1)
        else:
            Phe_dict = Phe[t].to_dict()
            if args.figGenes:
                fGenes =  Common.read_file(args.figGenes,mode="list",vals=[0])
                for gene in fGenes:
                    if gene not in Genes:
                        continue
                    single_lineRegress(Phe_dict,Exp.loc[gene].to_dict(),covar_df,outplot=True,gname=gene,tname=out)
            else:
                #run_lineRegress(Phe_dict,t,Exp,cut,f"{args.out}.{t}",args.loci)
                tasks.append((run_lineRegress,(Phe_dict,t,Exp,covar_df,cut,f"{args.out}.{t}",args.loci),{}))

    Common.run_parallel2(tasks, max_workers=args.threads, mode="func")
  


def read_Cov(cov_file):

    covar_df = None
    if cov_file is not None:
        covar_df = pd.read_csv(cov_file, sep="\t", header=None)
        n_pc = covar_df.shape[1] - 1
        covar_df.columns = ["sample"] + [f"PC{i+1}" for i in range(n_pc)]
        covar_df = covar_df.set_index("sample")

    return covar_df
    # --------- 样本对齐 ---------
    #samples = set(Phe.index) & set(gene_score_df.index)
    #if covar_df is not None:
    #    samples &= set(covar_df.index)
    #samples = sorted(samples)

    #if len(samples) == 0:
    #    raise ValueError("No common samples found among phenotype, gene_score, and covariates!")


def run_lineRegress(trait,trait_name,Exp,covar_df,cut,outfile,loci):
    # trait      : phe dict of single trait 
    # trait_name : trait name 
    # Exp        : express df of all genes


    print(f"## Run Association ... ... {trait_name}\n") 
    
    ### read Pos 
    Pos={}   
    if loci:
        of2=open(outfile+".rdata","w")
        of2.write("SNP\tChromosome\tPosition\tPvalue\n")

        for l in open(loci,"r"):
            ls=l.strip().split()
            mid = int((int(ls[1])+int(ls[2]))/2)
            Pos[ls[4]] = [ls[0],mid]


    of1=open(outfile+".txt","w")
    of1.write("Gene\tTrait\tlm_R2\tlm_b\tlm_p\tR\n")
        
    of3=open(outfile+".signal.txt","w")
    of3.write("Gene\tTrait\tlm_R2\tlm_b\tlm_p\tR\n")
    
    for gene in Exp.index.tolist():
        #model,pearson = single_lineRegress(trait,Exp.loc[gene].to_dict())
        model = single_lineRegress(trait,Exp.loc[gene].to_dict(),covar_df)
        r_squared = model.rsquared
        b = model.params[1]
        p_value = model.pvalues[1] if not str(model.pvalues[1]) == "nan" else  1.0
        R= r_squared**0.5
        #of1.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(gene,trait_name,r_squared,b,p_value,pearson[0],pearson[1]))
        of1.write(f"{gene}\t{trait_name}\t{r_squared}\t{b}\t{p_value}\t{R}\n")

        if float(p_value) <= cut:
            of3.write(f"{gene}\t{trait_name}\t{r_squared}\t{b}\t{p_value}\t{R}\n")

        if loci:
            of2.write(f"{gene}\t{Pos[gene][0]}\t{Pos[gene][1]}\t{p_value}\n")


    of1.close()
    of3.close()
    if loci:
        of2.close()


def single_lineRegress(Phe,Express,covar_df,outplot=False,gname=None,tname=None):

    rows = []
    for sam in Phe:
        if Phe[sam] == "NA":
            continue
        if sam not in Express:
            continue
        if covar_df is not None and sam not in covar_df.index:
            continue

        row = {
            "y": float(Phe[sam]),
            "Express": float(Express[sam])
        }

        if covar_df is not None:
            for c in covar_df.columns:
                row[c] = covar_df.loc[sam, c]

        #print(row)
        rows.append(row)

    df = pd.DataFrame(rows)
              
    #print(df)
    #xx = sm.add_constant(x)
    #model = sm.OLS(y,xx).fit()
    y = df["y"]
    X = df.drop(columns=["y"])
    X = sm.add_constant(X)

    model = sm.OLS(y, X).fit()
    #pearson = pearsonr(x, y)

    if not outplot:
        return model#,pearson

    r_squared = model.rsquared 
    p_value = model.pvalues[1] 


    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(x, y)
    ax.set_xlabel(gname, fontsize=12)
    ax.set_ylabel(tname, fontsize=12)
    ax.tick_params(axis='both', which='major', labelsize=10)

    # 画拟合线
    ax.plot(x, model.fittedvalues, color='red')

    # 在图上标注 R² 和 P
    ax.text(
        0.05, 0.95,
        f"R² = {r_squared:.3f}\nP = {p_value:.3e}",
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment='top',
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.7)
    )
    plt.savefig('TransAsso.%s_vs_%s.png' %(tname,gname),dpi=300, bbox_inches='tight')




if __name__ == "__main__":
    main()
