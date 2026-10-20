#!/usr/bin/env python3
import os
import sys
import Common
import Visualize
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import argparse

def get_args(args):
    parser = argparse.ArgumentParser(
    formatter_class = Common.CustomFormatter,
    description = '''
RunCmd :
    Funciton:  Performing Bulk-sequencing Analysis using Euclidean-distance or SNP-index.
    Writer:    Meng Minghui 
'''
    )
    parser.add_argument('-v','--vcf',help='vcf file',required=True,type=str)
    parser.add_argument('-t','--type',help='Use Euclidean-distance or SNP-index',type=str,default="ed",choices=["ed","index"],required=True)
    parser.add_argument('-o','--out',help='pre-name of out files',type=str,default="BSA")
    parser.add_argument('-s','--sign',help="Significance cut line: B = 1/n || F = 0.05/n || top1 = top 0.01 || top5 = top 0.05 || a self-defined float",nargs="+",default=["top1"])
    parser.add_argument('-step',help='Step-size of sliding window, kb',type=int,default=None)
    parser.add_argument('-wind',help='Wind-size of sliding window , kb',type=int,default=None)
    parser.add_argument('-fit',help='Fit method of sliding window',choices=["mean","linear"],type=str,default="mean")
    parser.add_argument('-minDepth',help='min Depht of each bulks',type=int,default=10)
    parser.add_argument('-minSite',help='min site of a window',type=int,default=10)
    parser.add_argument('-s1',help='Name of Segregation Bulk1',type=str,default=None,required=True)
    parser.add_argument('-s2',help='Name of Segregation Bulk2',type=str,default=None,required=True)
    parser.add_argument('-p1',help='Name of Parent Bulk1 ',type=str,default=None)
    parser.add_argument('-p2',help='Name of Parent Bulk2',type=str,default=None)
    
    parsed_args = parser.parse_args(args)
    print (parsed_args)
    return parsed_args 

def main(args=None):
    args=get_args(args)
    main_BSA(args)


def main_BSA(args):

    site_out= f"{args.out}.site.txt"

    if not os.path.exists(site_out):
        out1=open(site_out,"w")
        if args.type == "ed":
            out1.write("Chrom\tPos\tID\tRef\tAlt\ted2\ted4\tBulk1\tBulk2\n")
        elif args.type == "index":
            out1.write("Chrom\tPos\tID\tRef\tAlt\tGT_p1\tGT_p2\tabs_delta_index\tOrder\tdelta_index\tIndex1\tIndex2\tBulk1\tBulk2\n")
        name_index = {}

        print("# Start cal site BSA.... ")
        log = open(f"{args.out}.log","w")
        for line in open(args.vcf):
            if line.startswith("##"):
                continue
            elif line.startswith("#CHR"):
                for i, name in enumerate(line.strip().split()):
                    name_index[name] = i  
            else:
                info = line.strip().split()
                if args.type == "ed":
                    result = calED(info,name_index[args.s1],name_index[args.s2],args.minDepth)
                elif  args.type == "index":
                    result = calIndex(info,name_index[args.s1],name_index[args.s2],name_index[args.p1],name_index[args.p2],args.minDepth)
    
                if isinstance(result, str):
                    log.write(result+"\n")
                    continue
                else:
                    out1.write("\t".join(map(str,result)) + "\n") 
                #if result == None or len(result) == 0 :
                #    continue
                #out1.write("\t".join(map(str,result)) + "\n") 

    else:
        print("# result file {site_out} exist, remove it if you want re-run site-BSA")

    Site = pd.read_table(f"{args.out}.site.txt", sep="\t")
        
    if args.type == "ed":
        df_plot = Site.rename(columns={"ed2":"Value"})
        #print(df_plot)
        cutoff=Common.judge_significant(args.sign,Site,"ed2")

        Visualize.manhatan_fig_v1(df_plot,cutoff,"ed2",f"{args.out}.site")
        output_candidate(Site,"ed2",f"{args.out}.site",args.sign,cutoff)

        #winds
        if args.wind:
            Wind=sliding_window_fit(Site,f"{args.out}.wind{args.wind}k.txt",args.step*1000,args.wind*1000,args.minSite,method=args.fit,chrom_col='Chrom', pos_col='Pos', value_col='ed2')

            cutoff1=Common.judge_significant(args.sign,Wind,"Value")
            draw_combined(Site,Wind.fillna(0),cutoff1,f"{args.out}.wind{args.wind}k",posC="ed2")
            candidate=output_candidate(Wind,"Value",f"{args.out}.wind{args.wind}k",args.sign,cutoff1,mergeOut=True)


    elif args.type == "index":
        df_plot = Site.rename(columns={"abs_delta_index":"Value"})
        cutoff=Common.judge_significant(args.sign,Site,"abs_delta_index")

        Visualize.manhatan_fig_v1(df_plot,cutoff,"delta-index",f"{args.out}.site")
        output_candidate(Site,"abs_delta_index",f"{args.out}.site",args.sign,cutoff)

        if args.wind:
            Wind=sliding_window_fit(Site,f"{args.out}.wind{args.wind}k.txt",args.step*1000,args.wind*1000,args.minSite,method=args.fit,chrom_col='Chrom', pos_col='Pos', value_col='delta_index')

            Wind["Value_abs"] = Wind["Value"].abs()
            cutoff1=Common.judge_significant(args.sign,Wind,"Value_abs")
            cutoff2=cutoff1.copy()
            for s in cutoff1:
                cutoff2.append(-1*s)

            df_plot2 = Site.rename(columns={"delta_index":"Value"})
            draw_combined(df_plot2,Wind.fillna(0),cutoff2,f"{args.out}.wind{args.wind}k")
            candidate=output_candidate(Wind,"Value_abs",f"{args.out}.wind{args.wind}k",args.sign,cutoff1,mergeOut=True)
   


def output_candidate(df,colum,out,sign,cuts,mergeOut=False):
    for i,c in enumerate(cuts):
        df_c = df[df[colum] > cuts[i]]
        df_c.to_csv(f"{out}.candidate.{sign[i]}.txt",index=False,sep="\t")
        print(f"# Output candidate sites/regions to {out}.candidate.{sign[i]}.txt ")
        if mergeOut:
            merged_candidate=Common.merge_overlapping_intervals(df_c)
            merged_candidate.to_csv(f"{out}.candidate.{sign[i]}.merged.txt",index=False,sep="\t")


def draw_combined(Site,Wind,cuths,outfile,posC="Value",lineC="Value"):
    xmax=Site['Pos'].max()
    ymax=int(Site[posC].max())+1

    bins=list(Wind['Chrom'].unique())
    fig,axes = plt.subplots(len(bins),1,sharex=True,figsize=(2*len(bins),len(bins)))
    for i in range(len(bins)):
        name = bins[i]
        axes[i].text(xmax*1.05,ymax/2,name,fontsize=15,verticalalignment="top",horizontalalignment="left")
        axes[i].scatter(Site[Site['Chrom']==name]['Pos'],Site[Site['Chrom']==name][posC],s=0.5,c='lightgray')
        axes[i].plot(Wind[Wind['Chrom']==name]['Start'],Wind[Wind['Chrom']==name][lineC],c="firebrick")
        axes[i].spines['top'].set_color("white")
        axes[i].spines['right'].set_color("white")
        axes[i].xaxis.set_major_formatter(plt.FuncFormatter(millions))

        for c in cuths:
            axes[i].axhline(y=c,ls="--",c="navy",lw=1)

    plt.savefig(fname="%s.png" % (outfile),dpi=300)
    plt.show()



def sliding_window_fit(df,out,step_size,wind_size,minSite,method="mean",chrom_col='Chrom', pos_col='Pos', value_col='Value'):

    results = []
    out=open(out,"w")
    out.write(f"Chrom\tStart\tEnd\tValue\tNumber\n")
    total_winds=0
    vaild_winds=0

    for chrom in df[chrom_col].unique():
        chrom_data = df[df[chrom_col] == chrom].copy()
        chrom_data = chrom_data.sort_values(pos_col).reset_index(drop=True)
        
        if len(chrom_data) == 0:
            continue
            
        positions = chrom_data[pos_col].values
        values = chrom_data[value_col].values
        
        min_pos = positions.min()
        max_pos = positions.max()
        
        left_ptr = 0  
        window_id = 0
        window_start = 0
        
        while window_start <= max_pos:
            total_winds += 1
            window_end = window_start + wind_size
            
            while left_ptr < len(positions) and positions[left_ptr] < window_start:
                left_ptr += 1
            
            right_ptr = left_ptr
            while right_ptr < len(positions) and positions[right_ptr] < window_end:
                right_ptr += 1
            
            snp_indices = slice(left_ptr, right_ptr)
            snp_number = right_ptr - left_ptr
            
            fitted_value = None
            if snp_number >= minSite:
                vaild_winds += 1 
                fitted_value = fit_values(positions[snp_indices],values[snp_indices],method)
                    
            out.write(f"{chrom}\t{window_start}\t{window_end}\t{fitted_value}\t{snp_number}\n")
            #print(chrom,window_start,window_end,fitted_value,snp_number)
            results.append({
                'Chrom': chrom,
                'Start': window_start,
                'End': window_end,
                'Value': fitted_value,
                'Number': snp_number
            })

            
            window_start += step_size
            window_id += 1
    print(f"# Get {vaild_winds} vailded windows out of {total_winds}")

    return pd.DataFrame(results)
    

def fit_values(window_positions,window_values,method):
    if method == "mean":
        fitted_value = np.mean(window_values)
    
    elif method == "linear":
        if snp_number > 1:
            slope, intercept, _, _, _ = stats.linregress(window_positions, window_values)
            mid_pos = (window_start + window_end) / 2
            fitted_value = intercept + slope * mid_pos
        else:
            fitted_value = window_values[0]
    
    else:
        raise ValueError("method must be 'mean' or 'linear'")

    return fitted_value



def calED(info,s1,s2,minDepth):
    ds1 = info[s1].split(":")[1]
    ds2 = info[s2].split(":")[1]

    a1,a2 = map(int,ds1.split(","))
    b1,b2 = map(int,ds2.split(","))
    total =  a1 + a2 + b1 + b2
    if a1+a2 < minDepth or b1+b2 < minDepth:
        return f"Filter\t{info[0]}:{info[1]}\t#The Depth is less than the given threshold."
    else:
        ed2 = (float(a1)/(a1+a2)-float(b1)/(b1+b2))**2 + (float(a2)/(a1+a2)-float(b2)/(b1+b2))**2
        ed4 = ed2**2
        res= [info[0],info[1],info[2],info[3],info[4],ed2,ed4,ds1,ds2]
        return res


def calIndex(info, s1, s2, p1, p2, minDepth):
    """
    Calculate SNP-index and ΔSNP-index for BSA.
    
    Parameters:
        info: list of sample strings from VCF (GT:AD:DP:GQ:PL format)
        s1, s2: indices of bulk1 and bulk2
        p1, p2: indices of parent1 and parent2
        minDepth: minimum depth threshold for bulks
    
    Returns:
        (snp_index_1, snp_index_2, delta_index)
        Return None if any condition not satisfied.
    """

    # Parse parents
    p1_data = parse_sample(info[p1])
    p2_data = parse_sample(info[p2])
    
    if p1_data is None or p2_data is None:
        return f"Filter\t{info[0]}:{info[1]}\t#Invalid GT or AD information."
    
    gt_p1, ref_p1, alt_p1, dp_p1 = p1_data
    gt_p2, ref_p2, alt_p2, dp_p2 = p2_data

    # Normalize genotype separator
    gt_p1 = gt_p1.replace('|', '/')
    gt_p2 = gt_p2.replace('|', '/')

    # Parents must be homozygous
    if gt_p1 not in ("0/0", "1/1"):
        return f"Filter\t{info[0]}:{info[1]}\t#The genotype of parent 1 is not homozygous."
    if gt_p2 not in ("0/0", "1/1"):
        return f"Filter\t{info[0]}:{info[1]}\t#The genotype of parent 2 is not homozygous."

    # Parents must be different
    if gt_p1 == gt_p2:
        return f"Filter\t{info[0]}:{info[1]}\t#The genotypes of parent 1 and parent 2 are the same."

    # Parse bulks
    s1_data = parse_sample(info[s1])
    s2_data = parse_sample(info[s2])

    if s1_data is None or s2_data is None:
        return f"Filter\t{info[0]}:{info[1]}\t#Invalid GT or AD information."

    gt_s1, ref_s1, alt_s1, dp_s1 = s1_data
    gt_s2, ref_s2, alt_s2, dp_s2 = s2_data

    # Bulk depth filter
    if dp_s1 <= minDepth or dp_s2 <= minDepth:
        return f"Filter\t{info[0]}:{info[1]}\t#The Depth is less than the given threshold."

    # Determine which allele to track (parent1 allele)
    if gt_p1 == "0/0":
        snp_index_1 = ref_s1 / dp_s1
        snp_index_2 = ref_s2 / dp_s2
    else:
        snp_index_1 = alt_s1 / dp_s1
        snp_index_2 = alt_s2 / dp_s2

    delta_index = snp_index_1 - snp_index_2

    #out1.write("Chrom\tPos\tID\tRef\tAlt\tValue\tIndex1\tIndex2\tBulk1\tBulk2\n")
    order = "+" if delta_index > 0 else "-"
    res= [info[0],info[1],info[2],info[3],info[4],gt_p1,gt_p2,abs(delta_index),order,delta_index,snp_index_1,snp_index_2,f"{ref_s1};{alt_s1}",f"{ref_s2};{alt_s2}"]
    return res


def parse_sample(sample_str):
    fields = sample_str.split(':')
    if len(fields) < 2:
        return None
    
    gt = fields[0]
    ad = fields[1]
    # Disallow missing values
    if '.' in gt or '.' in ad:
        return None
    
    
    # Parse AD (must be biallelic: ref,alt)
    ad_parts = ad.split(',')
    if len(ad_parts) != 2:
        return None
    
    try:
        ref = int(ad_parts[0])
        alt = int(ad_parts[1])
        dp  = ref + alt
    except:
        return None
    
    return gt, ref, alt, dp


def millions(x,pos):
    return '{:1.1f}M'.format(x*1e-6)


if __name__ == '__main__':
    main()
