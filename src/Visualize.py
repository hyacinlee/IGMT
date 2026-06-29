#!/usr/bin/env python3 
import sys 
import Common
import argparse
import itertools
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt



def get_args(args):
    parser = argparse.ArgumentParser(description="Processing muti-traits GWAS with Emmax form vcf file",
                                    formatter_class=Common.CustomFormatter
        )

    base_group = parser.add_argument_group('Basic optional arguments')
    base_group.add_argument("-m","--model",required=True,help="Plot model: manha, pvg, stat",choices=["manha", "pvg", "stat"],default=None)
    base_group.add_argument("-o","--out",help="Prefix of outfiles",type=str,default="out")
    base_group.add_argument("-p","--phe",help="Phenotype file",type=str,default=None)
    base_group.add_argument("-f","--figsize",help="Phenotype file",nargs="+",default=[16,5])
    base_group.add_argument("-d","--dpi",help="DPI for png outfiles",type=int,default=300)
    base_group.add_argument("-c","--colors",help="colors in fig, one color per line in color file",type=str)
    base_group.add_argument("--ftype",help="output figure file types",type=str,choices=["png", "pdf"],default="png")

    manha_group = parser.add_argument_group('Manhantun optional arguments (manha)')
    manha_group.add_argument("--site_value",help="Gwas or Other site-value result file")
    manha_group.add_argument("--chrom",help="Chromosom name column",type=str,default="Chrom") 
    manha_group.add_argument("--pos",help="Maker pos column",type=str,default="Pos")
    manha_group.add_argument("--value",help="Value column",type=str,default="Value")
    manha_group.add_argument("--name",help="SNP name column",type=str,default="ID")  
    manha_group.add_argument("--sign",help="Significance cut line: B = 1/n || F = 0.05/n || top1 = top 0.01 || top5 = top 0.05  || a self-defined float",nargs="+",default=["B"])
    manha_group.add_argument("--log",help="change to log10 value",action="store_true",default=False)
    manha_group.add_argument("--sub",help="Draw only one chr",type=str,default=None)
    manha_group.add_argument("--local",help="Draw only local region be like chr:start-end",type=str,default=None)
    manha_group.add_argument("--hightlight",help="Hightlight maker ID in files, ids must be same with input",type=str,default=None)
    
    gvp_group = parser.add_argument_group('Correlation between genotype and phenotype optional arguments (pvg)')
    gvp_group.add_argument("--snps",help="Hightlight maker ID in files, ids must be same with input",type=str,default=None)
    gvp_group.add_argument("--gt",help="Traits assosciated GT file( *BasicGT.lst )",type=str,default=None)
    gvp_group.add_argument("--gt_start",help="The start col of GT file",type=int,default=6)
    gvp_group.add_argument("--gt_index",help="The col of maker ID name",type=int,default=1)
    gvp_group.add_argument("--group",help="sample in different group",type=str,default=None)

    stat_group = parser.add_argument_group('Visualize Basic stats using DataFrame file (stat)')
    stat_group.add_argument("--df",help="Input DataFrame file: header tsv file ",type=str,default=None)
    stat_group.add_argument("--plot",help="Choose to calculate by column or row",type=str,choices=["Col", "Row"],default="Col")
    stat_group.add_argument("--cor",help="Draw correlation heatmap (depending on --plot)",action="store_true",default=False)
    stat_group.add_argument("--pca",help="Perform and plot PCA (principal component analysis)",action="store_true",default=False)
    stat_group.add_argument("--box",help="Draw boxplot for data distribution using --xData and --yData",action="store_true",default=False)
    stat_group.add_argument("--lm_dot",help="Draw scatter plot using --xData and --yData columns",action="store_true",default=False)
    stat_group.add_argument("--exc",help="The non-numeric columns that needs to be excluded when run  --pca or --cor",nargs="*",default=None)
    stat_group.add_argument("--group_col",help="Groups col name for --pca or --dot",type=str,default=None)
    stat_group.add_argument("--group_color",help="Color  for groups in --box, --pca or  --dot",type=str,default=None)
    stat_group.add_argument("--group_min_rate",help="Filter groups less than group_min*Total number, also effective in -m pvg",type=float,default=0.01)
    stat_group.add_argument("--ttest",help="Chose t-test compare pairs like: g1,g2 g2,g3 g3,g4 ... ... when group > 3 ",nargs="*",default=None)
    stat_group.add_argument("--xData",help="Column name used as X-axis data",type=str,default=None)
    stat_group.add_argument("--yData",help="Column name used as Y-axis data",type=str,default=None)
    stat_group.add_argument("--xlab",help="Custom label for X-axis in plots",type=str,default=None)
    stat_group.add_argument("--ylab",help="Custom label for Y-axis in plots",type=str,default=None)

    parsed_args = parser.parse_args(args)
    print (parsed_args)
    return parsed_args
 

def main(args=None):

    args=get_args(args)

    figsize=(int(args.figsize[0]),int(args.figsize[1]))

    if args.model == "pvg":
        PvgMain(args)

    if args.model == "manha":
        ManhantanMain(args)

    if args.model == "stat":

        df=pd.read_csv(args.df, sep="\t")

        if args.plot == "Row":
            df = df.T

        if not args.xlab: 
            args.xlab = args.xData
        if not args.ylab:
            args.ylab=args.yData

        if args.group_color:
            colors = Common.read_file(args.group_color,mode="dict",vals=[1],keys=[0])
        else:
            colors = None

        if args.box:
            plot_boxplot_ttest(df,x=args.xData,y=args.yData,out=f"{args.out}.boxplot.{args.ftype}",
                               ttest_pairs=args.ttest,min_count=len(df)*args.group_min_rate,
                               figsize=figsize,dpi=args.dpi,xlab=args.xlab, ylab=args.ylab,colors=colors)
        if args.cor:
            plot_cor(args.df,f"{args.out}.cor_{args.plot}.{args.ftype}",f"{args.out}.cor_matrix.tsv",cor_data=args.plot,excluded=args.exc,figsize=figsize,dpi=args.dpi)

        if args.lm_dot:
            plot_lm_dot(df,x=args.xData,y=args.yData,out=f"{args.out}.cor_{args.plot}.{args.ftype}",
                       xlab=args.xlab, ylab=args.ylab,figsize=figsize,dpi=args.dpi)



def plot_lm_dot(df,x,y,out,xlab,ylab,figsize=(8,5),dpi=300):
    #print(df)
    x=df[x].tolist()
    y=df[y].tolist()
             
    xx = sm.add_constant(x)
    model = sm.OLS(y,xx).fit()

    r_squared = model.rsquared 
    p_value = model.pvalues[1] 

    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(x, y)
    ax.set_xlabel(xlab, fontsize=10)
    ax.set_ylabel(ylab, fontsize=10)
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
    plt.savefig(out, dpi=dpi, bbox_inches='tight')



def plot_cor(infile,outfig,outtsv,cor_data="Col",excluded=None,figsize=(10,10),dpi=300):

    df = pd.read_csv(infile, sep="\t", index_col=0)
    if excluded:
        df = df.drop(columns=excluded, errors="ignore")

    if cor_data == "Col":
        correlation_matrix = df.corr()
    elif cor_data == "Row":
        correlation_matrix = df.T.corr()
    else:
        raise ValueError("--plot must be 'Col' or 'Row'")


    if correlation_matrix.shape[0] < 2:
        raise ValueError("Not enough data for clustering — check that your matrix has multiple genes/samples.")
    
    g = sns.clustermap(correlation_matrix, cmap="coolwarm", figsize=figsize)
    plt.title("Correlation clustering plot")

    plt.savefig(outfig, dpi=dpi, bbox_inches="tight")
    correlation_matrix.to_csv(outtsv,sep="\t")





def ManhantanMain(args):
    
    name_col= {args.chrom:"Chrom",args.pos:"Pos",args.value:"Value",args.name:"ID"}
    df = pd.read_table(args.site_value, sep="\t")
    df = df.rename(columns=name_col)
    df = df.dropna(how="any", axis=0)
  
    df,cut,yname = clean_logs_cut(df,args.log,args.sign,args.value)

    colors=["#3274a1", "#e1812c"] if not args.colors else Common.read_file(args.colors,mode="list",vals=[0],skipAnnot=False)
    print(f"# Manhantun chromosomes color: {colors}")
    hightlight = None if not args.hightlight else Common.read_file(args.hightlight,mode="dict",keys=[0],vals=[1,2])
    print(f"# Hightlight SNP markers: {hightlight} ")

    #if args.local:
    #    local_plot(df,cut,yname,args.out,args.local,hightlight,colors)

    manhatan_fig_v1(df,cut,yname,args.out,args.sub,args.local,hightlight,colors,args.figsize,args.dpi,args.ftype)

 
def manhatan_fig_v1(df, cut, ylab, out, sub=None,local=None,hightlight=None,colors=["#3274a1", "#e1812c"],figsize=[16,5],dpi=300,ftype="png"):
    """
    绘制曼哈顿图，支持高亮SNP注释、染色体子集和多色，支持双向绘制（根据Order列）。
    df: 包含 Chrom, Pos, Value, ID 四列的数据框（如果有Order列则用于双向绘图）
    cut: 显著性阈值线，可以是一个或多个值的数组（会自动镜像到负方向）
    ylab: y轴标签
    out: 输出文件名（不带扩展名）
    sub: 若指定，只画某条染色体
    hightlight: 字典，key为SNP ID，value为注释文本
    colors: 颜色数组用于交替染色
    figsize:  长宽比
    dpi: 输出图像分辨率
    ftype：输出类型
    """
    df = df.copy()
    
    has_direction = 'Order' in df.columns
    
    y_max = df['Value'].max()

    if sub:
        df = df[df['Chrom'] == sub]

    if local:
        #print(df)
        subchr,se=local.split(":")
        start=int(se.split("-")[0])
        end  =int(se.split("-")[1])
        print(f"# Draw plot map in {subchr} from {start} to {end}") 
        df =df[(df['Chrom'] == subchr) & (df['Pos'] <= end) & (df['Pos'] >= start)]
        out = f"{subchr}_{start}_{end}"

    print(df)

    df['Chrom'] = df['Chrom'].astype(str)
    chrom_order = Common.get_sorted_chromosomes(df['Chrom'].unique())
    df['Chrom'] = pd.Categorical(df['Chrom'], categories=chrom_order, ordered=True)
    df = df.sort_values(['Chrom', 'Pos'])

    if has_direction:
        df['Value'] = df.apply(lambda row: -row['Value'] if row['Order'] == '-' else row['Value'], axis=1)


    df['ind'] = range(len(df))
    df_grouped = df.groupby('Chrom')

    fw=int(figsize[0])
    fh=int(figsize[1])
    print(fw,fh)
    fig, ax = plt.subplots(figsize=(fw,fh))

    x_labels = []
    x_labels_pos = []
    
    for i, (name, group) in enumerate(df_grouped):
        print(f"# Plotting Chrom: {name} into manhattan figure")
        ax.scatter(group['ind'], group['Value'],
                   color=colors[i % len(colors)], s=10)
        x_labels.append(name)
        x_labels_pos.append((group['ind'].iloc[0] + group['ind'].iloc[-1]) / 2)


    for c in cut:
        ax.axhline(y=c, color='gray', linestyle='--', linewidth=1)
        if has_direction:
            ax.axhline(y=-c, color='gray', linestyle='--', linewidth=1)

    if hightlight:
        highlight_ids = list(hightlight.keys())
        highlight_df = df[df['ID'].isin(highlight_ids)]
    
        for _, row in highlight_df.iterrows():
            snp_id = row['ID']
            x = row['ind']
            y = row['Value']
            label, color = hightlight.get(snp_id, [snp_id, 'red'])
            
            ax.scatter(x, y, color=color, s=10, zorder=10)
            

            if label != "-":
                y_text = y + (0.8 if y >= 0 else -0.8)  
                ax.annotate(label,
                           xy=(x, y), xycoords='data',
                           xytext=(x, y_text), textcoords='data',
                           arrowprops=dict(arrowstyle='-', color=color, lw=0.8),
                           ha='center', va=('bottom' if y >= 0 else 'top'), 
                           fontsize=8, color=color, rotation=45)

    # 设置轴
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel(ylab)
    ax.set_title(out)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


    # 如果有方向信息，调整y轴范围并添加零线
    if has_direction:
        y_max = max(df['Value'].max(), abs(df['Value'].min()))
        ax.set_ylim(-y_max*1.1, y_max*1.1)
        ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    else:
        ax.set_ylim(0, y_max*1.1)

    ticks = ax.get_yticks()
    ax.set_yticklabels([f"{abs(tick):g}" for tick in ticks])
    #ax.set_yticks(ticks)  # 第一步：设置位置
    #if labels is None:
        #labels = [f"{abs(tick):g}" for tick in ticks]
        #ax.set_yticklabels(labels)  # 第二步：设置标签


    plt.tight_layout()
    if ftype == "png":
        plt.savefig(f"{out}.manhattan.png", dpi=dpi)
    else:
        plt.savefig(f"{out}.manhattan.pdf", dpi=dpi)

    print(f"# Finish output manhattan plot in {out}.manhattan")



def manhatan_fig_v0(df, cut, ylab, out, sub=None, hightlight=None,colors=["#3274a1", "#e1812c"],dpi=300):
    """
    绘制曼哈顿图，支持高亮SNP注释、染色体子集和多色。
    df: 包含 Chrom, Pos, Value, ID 四列的数据框
    cut: 显著性阈值线，可以是一个或多个值的数组
    ylab: y轴标签
    out: 输出文件名（不带扩展名）
    sub: 若指定，只画某条染色体
    hightlight: 字典，key为SNP ID，value为注释文本
    colors: 颜色数组用于交替染色
    """
    df = df.copy()
   
    # 如果指定了子集染色体
    if sub:
        df = df[df['Chrom'] == sub]

    # 排序染色体
    df['Chrom'] = df['Chrom'].astype(str)
    chrom_order = sorted(df['Chrom'].unique(), key=chrom_sort_key)
    df['Chrom'] = pd.Categorical(df['Chrom'], categories=chrom_order, ordered=True)
    df = df.sort_values(['Chrom', 'Pos'])

    # 计算位置偏移
    df['ind'] = range(len(df))
    df_grouped = df.groupby('Chrom')

    # 绘图
    fig, ax = plt.subplots(figsize=(16, 8))

    x_labels = []
    x_labels_pos = []
    
    for i, (name, group) in enumerate(df_grouped):
        print(f"# Plotting Chrom: {name} into manhandun figure")
        ax.scatter(group['ind'], group['Value'],
                   color=colors[i % len(colors)], s=10)
        x_labels.append(name)
        x_labels_pos.append((group['ind'].iloc[0] + group['ind'].iloc[-1]) / 2)

    # 添加cut线
    for c in cut:
        ax.axhline(y=c, color='gray', linestyle='--', linewidth=1)

    # 高亮 SNP 和注释
    if hightlight:
        highlight_ids = list(hightlight.keys())
        highlight_df = df[df['ID'].isin(highlight_ids)]
    
        # 单独处理每个高亮 SNP（因颜色不同）
        for _, row in highlight_df.iterrows():
            snp_id = row['ID']
            x = row['ind']
            y = row['Value']
            y_text = y + 0.8  # 文本高度
            label, color = hightlight.get(snp_id, [snp_id, 'red'])
    
            # 画点
            ax.scatter(x, y, color=color, s=25, zorder=10)
    
            # 画注释和线
            ax.annotate(label,
                        xy=(x, y), xycoords='data',
                        xytext=(x, y_text), textcoords='data',
                        arrowprops=dict(arrowstyle='-', color=color, lw=0.8),
                        ha='center', va='bottom', fontsize=8, color=color, rotation=45)

    # 设置轴
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel(ylab)
    ax.set_title(out)

    plt.tight_layout()
    plt.savefig(f"{out}.manhantun.png",dpi=dpi)
    print(f"# Finish ouptut manhantun plot in {out}.manhantun.png")


def clean_logs_cut(df,log,sign,name):
    cut = Common.judge_significant(sign,df)  
    
    #yname="-log10(P)"
    #cut = -np.log10(cut)

    #if not log:
    #    df["Value"]=1.0/10**df["Value"]
    #    yname=name
    #    cut = 1/10**float(cut)
    #else:
    #    df['Value'] = -np.log10(df['Value'])
    if log:
        df['Value'] = -np.log10(df['Value'])
        cut = -np.log10(cut)
        yname="-log10(P)"
    else:
        yname = name    


    return df,cut,yname

def chrom_sort_key(chrom):
    # 提取数字部分作为排序主键；X/Y/MT等排在最后
    try:
        return int(chrom)
    except:
        # 自定义顺序：X=23, Y=24, MT=25，其他未知为99
        chrom_map = {'X': 23, 'Y': 24, 'MT': 25}
        return chrom_map.get(chrom.upper(), 99)



def PvgMain(args):
    
    GTs,header = Common.read_file(args.gt,mode="dict",keys=[args.gt_index-1],vals=[],header=True)
    snps = Common.read_file(args.snps,mode="list",vals=[0])
    Phe = Common.read_file(args.phe,mode="dict",keys=[0],vals=[2])
   
    df_inf = pd.DataFrame.from_dict(Phe, orient='index', columns=['Phenotype'])
    df_inf.index.name="Sample"
    #df_inf = pd.to_numeric(df_inf, errors="coerce").dropna()

    
    if args.group:
        sample_group= Common.read_file(args.group,mode="dict",keys=[0],vals=[1])
        df_inf=Common.updata_df(df_inf,"Group",sample_group,missing="None")
        df_draw=df_inf[df_inf["Group"] != 'None']
        plot_boxplot_ttest(df_draw,'Group', 'Phenotype',"Phe_by_Group.pdf")

    print(df_inf)  

    for snp in snps:
        sample_gt = dict(zip(header[args.gt_start-1:],GTs[snp][args.gt_start-1:]))
        df_inf = Common.updata_df(df_inf,snp,sample_gt,missing="None")
        df_draw1=df_inf[df_inf[snp] != 'NN']
        df_draw1=df_inf[df_inf['Phenotype'] != 'NA']
        #min_count = len(df_draw1)*0.05
        plot_boxplot_ttest(df_draw1,snp, 'Phenotype',f"{snp}.Phe_by_GT.pdf",min_count=len(df_draw1)*args.group_min_rate)
        if args.group:
            df_draw2=df_draw1[df_draw1['Group'] != 'None']
            plot_stacked_bar(df_draw2,f"{snp}.GT_by_group.pdf",value=snp,min_count=len(df_draw2)*args.group_min_rate)

    df_inf.to_csv(f"{args.out}.phe.gt.info.xls",sep="\t",index=True)


def plot_stacked_bar(df,out,value="Phenotype"):
    """
    绘制按照 Group 划分的 Value 堆叠柱状图。
    
    参数:
    df : pd.DataFrame
        包含三列：Sample、Group、Value
    """
    # 确保列名一致
    if not set(['Group', value]).issubset(df.columns):
        raise ValueError(f"DataFrame 必须包含 'Group' 和 {value}三列")

    # 将 Value 转换为计数：每个 Group 内的 Value 出现次数
    count_df = df.groupby(['Group', value]).size().unstack(fill_value=0)

    # 绘图
    ax = count_df.plot(kind='bar', stacked=True, figsize=(10, 6), colormap='tab20')

    plt.title(f'{value} by Group')
    plt.xlabel('Group')
    plt.ylabel('Count')
    plt.legend(title=value, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(out, format='pdf')


def plot_boxplot_ttest(df,x,y,out,ttest_pairs=None,figsize=(5, 5),dpi=300,xlab=None,ylab=None,colors=None,min_count=1):

    # === 数据检查 ===
    assert x in df.columns and y in df.columns, f"DataFrame must contain '{x}' and '{y}' columns"
    df = df.dropna(subset=[x, y]).copy()
    df[x] = df[x].astype(str)
    df[y] = pd.to_numeric(df[y], errors='coerce')

    # === 统计每组样本数，并按数量排序 ===
    group_counts = df.groupby(x)[y].count().sort_values(ascending=False)
    valid_groups = group_counts[group_counts >= min_count].index.tolist()

    # 如果提供了 colors，则仅保留 colors 中有定义的组（但不重排颜色）
    if colors is not None:
        valid_groups = [g for g in valid_groups if g in colors]

    df = df[df[x].isin(valid_groups)]

    if len(valid_groups) == 0:
        raise ValueError(f"No valid groups left to plot (min_count={min_count})")

    print(f"# Groups kept (n >= {min_count}, sorted by count desc): {valid_groups}")
    print(f"# Groups removed (n < {min_count}): {[g for g in group_counts.index if g not in valid_groups]}")

    gt_list = valid_groups  # 按样本数从大到小排列

    # === 绘图 ===
    plt.figure(figsize=figsize)
    if colors is not None:
        # 用户自定义颜色映射，不改变
        palette = colors
        ax = sns.boxplot(x=x, y=y, data=df, width=0.4,order=gt_list, palette=palette)
    else:
        ax = sns.boxplot(x=x, y=y, data=df, width=0.4,order=gt_list)

    # === 修改 x 轴标签，添加 n 值 ===
    n_per_group = df.groupby(x)[y].count()
    new_labels = [f"{g}\n(n={n_per_group.get(g, 0)})" for g in gt_list]
    ax.set_xticklabels(new_labels)

    # === t-test ===
    if ttest_pairs:
        pairs = [p.split(",") for p in ttest_pairs]
    elif len(gt_list) > 1 and len(gt_list) <= 3:
        pairs = list(itertools.combinations(gt_list, 2))
    else:
        pairs = []

    if pairs:
        print(f"# Run t-test for pairs: {pairs}")
        ymax = df[y].max()
        h = (df[y].max() - df[y].min()) * 0.05
        curr_y = ymax + h
        x_pos = {gt: i for i, gt in enumerate(gt_list)}

        for (gt1, gt2) in pairs:
            if gt1 not in df[x].unique() or gt2 not in df[x].unique():
                print(f"Skip invalid pair: {gt1}, {gt2}")
                continue
            group1 = df[df[x] == gt1][y]
            group2 = df[df[x] == gt2][y]
            print(group1)
            print(group2)
            stat, p = ttest_ind(group1, group2, equal_var=False)
            print(stat,p)
            #stat, p = ttest_ind(group1, group2)
            label = get_sig_label(p)

            x1, x2 = x_pos[gt1], x_pos[gt2]
            ax.plot([x1, x1, x2, x2],
                    [curr_y, curr_y + h, curr_y + h, curr_y],
                    lw=1.5, c='k')
            ax.text((x1 + x2) / 2, curr_y + h + 0.01, label,
                    ha='center', va='bottom', fontsize=10)
            curr_y += h * 2

    plt.title(f"Boxplot of {y} by {x}")
    plt.xlabel(xlab)
    plt.ylabel(ylab)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(out, dpi=dpi)
    plt.close()
    print(f"Boxplot with t-test saved to {out}")


def get_sig_label(p):
    if p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    else:
        return 'ns'


if __name__ == '__main__':

    main()
