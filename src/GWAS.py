#!/usr/bin/env python3
import Common
import Visualize
import ParaRun
import pandas as pd
import numpy as np
import gzip
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import argparse

def get_args(args):
    parser = argparse.ArgumentParser(description="Processing muti-traits GWAS with Emmax form vcf file",
                                    formatter_class=Common.CustomFormatter
        )
    parser.add_argument("-v", "--vcf",required=True, help="Path of the input vcf file (.vcf or .vcf.gz)")
    parser.add_argument("-p", "--phe",required=True, help="Path of the input phenotype file")
    parser.add_argument("-c","--cov",help="Covariance, could be PCA, structure or other files",type=str) 
    parser.add_argument("-o","--out",help="Prefix of outfiles",type=str,default="out")
    parser.add_argument("-m","--maf",help="Exclude variants with minor allele frequency lower than this",type=float,default="0.05")
    parser.add_argument("-s","--sign",help="Significance cut line: B = 1/n || F = 0.05/n || a self-defined float",nargs="+",default=["B"])
    parser.add_argument("-t","--threads",help="The number of parallel threads to perform GWAS",type=int,default=12)
    parser.add_argument("--split",help="Run GWAS separately by chromosome, specify the chromosome list",type=str)  
    parser.add_argument("--emmax",help="Path of the emmax dir,which include emmax-kin-intel64 and emmax-intel64 ",default="/opt/software/EMMAX/")
    parser.add_argument("--plink",help="Path of the plink v1.9",default="/opt/software/plink1.9/plink")
     
    parsed_args = parser.parse_args(args)
    print (parsed_args)
    return parsed_args
  
def main(args=None):
    args=get_args(args)

    Common.check_path_exists(f"{args.emmax}/emmax-intel64")
    Common.check_path_exists(f"{args.emmax}/emmax-kin-intel64")
    Common.check_path_exists(f"{args.plink}")

    samples=process_vcf_phe(args)

    if args.cov:
        cov_list=process_covariance(args.cov,samples)

    traits=process_gwas(args)

    with open("GWAS.traits","w") as f:
        for t in traits:
            f.write(f"{t}\n")

    snps=Common.read_file("GWAS.snpID","dict",vals=[0,1],keys=[2])

    significant=Common.judge_significant(args.sign,snps)

    with open("GWAS.threshold","w") as f:
        f.write(f"{significant}\n")
 
    if not Path("GWAS.outputs.done").exists():
        process_outputs(traits,significant,args.threads,len(samples))

    if not Path("GWAS.manhantun.done").exists():
        process_manhantun(traits,significant)


def process_manhantun(traits,significant):
    tasks=[]
    for trait in traits:
        df = pd.read_table(f"Traits.{trait}.emmax.ps.tsv.gz", sep="\t",compression="gzip")
        df['Value'] = -np.log10(df['Value'])
        cut = -np.log10(significant)
        #tasks.append((Visualize.manhatan_fig_v1,(df,cut,"-log(P)",f"Traits.{trait}"),{}))
        Visualize.manhatan_fig_v1(df,cut,"-log(P)",f"Traits.{trait}")
        #Visualize.ManhantanMain(f"Traits.{trait}.emmax.ps.tsv",out=f"Traits.{trait}.emmax.ps.manhantun",sign=significant)
    #Common.run_parallel(tasks, max_workers=5, mode="func")
    Common.run_command(f"touch GWAS.manhantun.done")



def process_outputs(traits,significant,threads,sampel_number):
    tasks=[]
    for trait in traits:
        #process_trait_output(trait,snps,significant,sampel_number)
        tasks.append((process_trait_output,(trait,significant,sampel_number),{}))

    Common.run_parallel2(tasks, max_workers=4, mode="func")
    #Common.run_parallel_thread(tasks, max_workers=threads)

    Common.run_command(f"touch GWAS.outputs.done")

    #return significant


def process_trait_output(trait,significant,sampel_number):
    print(f"# Output {trait} result...")

    #snps = Common.read_file("GWAS.snpID","dict",vals=[0,1],keys=[2])

    snps=Common.read_file("GWAS.snpID","dict",vals=[0,1],keys=[2])

    #with open(f"Traits.{trait}.emmax.ps", "r") as f,open(f"Traits.{trait}.emmax.ps.tsv", "w") as f1,open(f"Traits.{trait}.emmax.ps.significant.tsv", "w") as f2:
    with open(f"Traits.{trait}.emmax.ps", "r") as f, gzip.open(f"Traits.{trait}.emmax.ps.tsv.gz", "wt") as f1,open(f"Traits.{trait}.emmax.ps.significant.tsv", "w") as f2:
        f1.write("Chrom\tPos\tID\tValue\tBeta\tSE\tR2\n")
        f2.write("Chrom\tPos\tID\tValue\tBeta\tSE\tR2\n")

        #data=[]
        for l in f:
            ls = l.strip().split("\t")
            #chrom, pos = snps.loc[ls[0], ["chrom", "pos"]]
            if float(ls[3]) == 0.0 : # Avoid invalid p-values , update 20260127
                ls[3] = 0.99
            if float(ls[3]) == 1.0 : 
                ls[3] = 0.99
            chrom, pos = snps[ls[0]]
            if float(ls[2]) ==0.0:
                t = 0
            else:
                t = float(ls[1])/float(ls[2])
            r2 = t**2/(t**2+sampel_number-2)
            
            f1.write(f"{chrom}\t{pos}\t{ls[0]}\t{ls[3]}\t{ls[1]}\t{ls[2]}\t{r2}\n")
            if float(ls[3]) <= significant[0]:
                f2.write(f"{chrom}\t{pos}\t{ls[0]}\t{ls[3]}\t{ls[1]}\t{ls[2]}\t{r2}\n")


def process_gwas(args):
    data,traits,samples = Common.read_matrix_data(f"GWAS.phe","h")

    print (f"# Run GWAS with population size: {len(samples)} with {len(traits)} Traits")
    print (f"# Traits for run: {traits}")

    cmds=[]
    chr_list=[]
    for trait in traits:

        trait_phes,samples = Common.sort_dict_by_list(data[trait],samples)
        Common.write_list_to_file(zip(samples,samples,trait_phes),f"Traits.{trait}.phe",vals=[])

        if args.split:
            chr_list=Common.read_file(args.split,mode="list",vals=[0],header=False)
            Common.run_command(f"mkdir -p Traits.{trait}.split")
            if args.cov:
                cs=[f"{args.emmax}/emmax-intel64 -v -d 10 -t GWAS.split/{c} -k GWAS.aBN.kinf -p Traits.{trait}.phe -o Traits.{trait}.split/{c}.emmax -c GWAS.covariance" for c in chr_list]
            else:
                cs=[f"{args.emmax}/emmax-intel64 -v -d 10 -t GWAS.split/{c} -k GWAS.aBN.kinf -p Traits.{trait}.phe -o Traits.{trait}.split/{c}.emmax" for c in chr_list]
            cmds+=cs
        else:
            if args.cov:
                cmds.append(f"{args.emmax}/emmax-intel64 -v -d 10 -t GWAS -k GWAS.aBN.kinf -p Traits.{trait}.phe -o Traits.{trait}.emmax  -c GWAS.covariance")
            else:
                cmds.append(f"{args.emmax}/emmax-intel64 -v -d 10 -t GWAS -k GWAS.aBN.kinf -p Traits.{trait}.phe -o Traits.{trait}.emmax")
    
    Common.write_list_to_file(cmds,"./cmd.GWAS.sh")
    ParaRun.runlocal("./cmd.GWAS.sh",lines=1,threads=args.threads)


    if not len(chr_list) == 0:
        cmds=[]
        for trait in traits:
            if not Common.check_path_exists(f"Traits.{trait}.emmax.ps",mode = "w"):
                result_list =  [f"Traits.{trait}.split/{c}.emmax.ps" for c in chr_list]
                cmds.append(f"cat {' '.join(result_list)} > Traits.{trait}.emmax.ps && rm -rf Traits.{trait}.split/*emmax.ps")
            else:
                print(f"GWAS Result file: Traits.{trait}.emmax.ps exists, skip...." )
        if len(cmds) > 0 :
            Common.run_parallel2(cmds,max_workers=args.threads, mode="cmd", capture_output=True)

    return traits


def process_vcf_phe(args):

    # out 
    vcf_samples=[]
    vcf_open = gzip.open if args.vcf.endswith(".gz") else open

    if not Path("GWAS.snpID").exists():
        with vcf_open(args.vcf,"rt") as f1, open("GWAS.snpID","w") as fo:
            for line in f1:
                ls = line.strip().split("\t")
                if line.startswith("#CHR"):
                    vcf_samples = line.strip().split("\t")[9:]
                elif not line.startswith("#"):
                    if not ls[2] == ".":
                        fo.write(f"{ls[0]}\t{ls[1]}\t{ls[2]}\n")
                    else:
                        print(f"# There is no SNP-ID infomation cols in {args.vcf}, you can use the Vcftools module to add it")
                        exit(1)
    else:
        print(f"# File GWAS.snpID exists, passing") 
        with vcf_open(args.vcf,"rt") as f1:
            for line in f1:
                ls = line.strip().split("\t")
                if line.startswith("#CHR"):
                    vcf_samples = line.strip().split("\t")[9:]
                    break


    (phe_data,headinfo)=Common.read_file(args.phe,mode="dict",keys=[0],header=True)
    (phe_data,samples)=Common.sort_dict_by_list(phe_data,vcf_samples)

    Common.write_list_to_file(phe_data,"GWAS.phe",vals=[],header=headinfo)

    if Path("GWAS.prepared.done").exists():
        old_samples=Common.read_file("GWAS.sample","list",vals=[0])
        if sorted(old_samples) == sorted(samples):
            print(f"# The old samples equl the news samples. passing GWAS preparing")
            Common.write_list_to_file(phe_data,"GWAS.phe",vals=[],header=headinfo)
            return samples
        else:
            print(f"# GWAS.prepared.done exist,but the samples of pre_run is not equl with this run, please remove GWAS.prepared.done and rerun this commond!")
            exit(1)
    else:
        Common.write_list_to_file(phe_data,"GWAS.sample",vals=[0,0])
        #Common.write_list_to_file(phe_data,"GWAS.phe",vals=[],header=headinfo)
        Common.run_command(f"{args.plink} --vcf {args.vcf} --allow-extra-chr --make-bed --out GWAS --keep GWAS.sample --maf {args.maf}")
        if args.split:
            chr_list=Common.read_file(args.split,mode="list",vals=[0],header=False)
            Common.run_command(f"mkdir -p  GWAS.split")
            cmds=[f"{args.plink} --bfile GWAS --out GWAS.split/{c} --allow-extra-chr --chr {c} --output-missing-genotype 0 --recode 12 transpose" for c in chr_list]
            Common.run_parallel2(cmds,max_workers=8, mode="cmd", capture_output=True)
            #print(result)

        Common.run_command(f"{args.plink} --bfile GWAS --out GWAS --allow-extra-chr --output-missing-genotype 0 --recode 12 transpose")
        Common.run_command(f"{args.emmax}/emmax-kin-intel64 -v -x -d 10 GWAS")   # Generate kinship matrix
        Common.run_command(f"touch GWAS.prepared.done")  
        return samples


def process_covariance(cov_file,samples):
    """Process structure/PCA file"""

    cov_data=Common.read_file(cov_file,mode="dict",keys=[0],vals=[]) # read cov without header
    #print(cov_data)
    cov_data_new=[] 
    for s in samples:
        ls = cov_data[s]
        ls.insert(1,1)
        ls.insert(1,ls[0])
        cov_data_new.append(ls)

    Common.write_list_to_file(cov_data_new,f"GWAS.covariance",vals=[]) 

if __name__ == "__main__":
    main()
