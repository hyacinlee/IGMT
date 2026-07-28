#!/usr/bin/env python3 
import re
import os
import sys 
import argparse
import subprocess
import logging
import numpy as np
import pandas as pd
from scipy import stats
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from concurrent.futures import ProcessPoolExecutor, as_completed


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawTextHelpFormatter):
    def _format_action_invocation(self, action):
        if not action.option_strings:
            return super()._format_action_invocation(action)
        return ', '.join(action.option_strings)

    def _get_help_string(self, action):
        help_str = super()._get_help_string(action)
        
        if action.required:
            help_str += ' (required)'

        #if action.default:
        #    help_str = help_str.replace("(default: None)", "")
        
        return help_str



# File - List & Dict read/write  Function 
def updata_df(df,index_name,data_dict,missing="None"):
    df[index_name] = df.index.map(data_dict)
    df[index_name].fillna(missing, inplace=True)
    #new_order = [c for c in df.columns if c != 'B'] + ['B']
    return df

def read_file(infile,mode="list",vals=[],keys=[],header=False,sep="\t",noSplit=False,skipAnnot=True):


    result = [] if mode == "list" else {}
    head_info=[]
    print(f"# Read file {infile} to {mode},vals={[vals]},keys={[keys]},skip={skipAnnot}")
    
    with open(infile,"r") as inf:
        for i, line in enumerate(inf):
            if skipAnnot and line.startswith("#"):    # start with #
                continue

            vv = return_vals(vals,line,sep,noSplit)

            if header and i==0:
                head_info=vv
                continue

            kv = return_vals(keys,line,sep,False)

            if mode == "list":
                result.append(vv)
            elif mode == "dict":
                key = kv
                if type(key) == list:
                    key = "###".join(kv)
                result[key] = vv

    if header==False:
        return result
    else:
        return (result,head_info)


def read_file_accumulateDict(infile,vals=[0],key1=[1],key2="no",sep="\t"):

    result={}
    print(f"# Read file {infile} to accumulateDict,vals={[vals]},key1={[key1]},key2={[key2]}")
    with open(infile,"r") as inf:
        for line in inf:
            if line.startswith("#"):    # start with #
                continue
            datas = line.strip().split(sep)
            vv = return_vals(vals,line,sep,False)
            #print(vv)

            #print(key1)
            key1vv = return_vals(key1,line,sep,False)
            #

            if not key2 =="no":
                key2vv = return_vals(key2,line,sep,False)
                #print(key2)
            else:
                key2vv="no"
            #print(key1vv,key2vv,vv)

            result=accumulateDict(result,vv,key1vv,key2vv)
            
        return result



def return_vals(vals,line,sep="\t",noSplit=False):
    if noSplit: # not split ,return whole line 
        return line.strip()
    else:
        datas = line.strip().split(sep)
        for i in vals:
            if i > len(datas):
                print(f"#Error: index {i} not in infoline:{line}")

        vv=""
        if len(vals)==1:            # only one, return whole line in str
            vv = datas[vals[0]]
        elif len(vals)==0:          # return all    cols in lsit 
            vv=datas
        else:                       # return parted cols in lsit 
            vv=[datas[i] for i in vals]
    return vv



def write_list_to_file(mylist,output,vals=None,header=False,sep="\t"):
    # write a list to file
    # vals = None : a simple list 
    # vals != None: a two dimom list

    with open(output, 'w') as f:

        if vals == None:
            if header:
                f.write(header+"\n")
            for s in mylist:
                f.write(s+"\n")
        else:
            if header:
                datas = header if len(vals) == 0 else [ header[i] for i in vals ]
                f.write(sep.join(datas)+"\n")
            for row in mylist:
                row = list_type(row,"str")
                datas = row if len(vals) == 0 else [ row[i] for i in vals ]
                f.write(sep.join(datas)+"\n")


def sort_dict_by_list(mydict,sort_list):
    # sort and filter a dict by a sorted list
    # return the sorted dict value and list 
    vs=[]
    ks=[]
    for s in sort_list:
        if s in mydict:
            ks.append(s)
            vs.append(mydict[s])
    return vs,ks


def read_matrix_data(infile,order):
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
                try:
                    v = ls[i]
                except IndexError:
                    print(f"i = {i}")
                    print(f"len(ls) = {len(ls)}")
                    print(ls)
                    print(l.rstrip())
                    raise


                if order == "h":  ##  H-line  as key
                    P=none_rep_twoDict(P,v,k1,k2) 
                if order == "v":  ##  V-line  as key
                    P=none_rep_twoDict(P,v,k2,k1)
    return P,Hlist,Vlist

def none_rep_twoDict(myDict,val,keyA,keyB):
    if keyA not in myDict:
        myDict[keyA]={}
    myDict[keyA][keyB]=val

    return myDict



def accumulateList(mylist,val):
    if val not in mylist:
        mylist.append(val)
    return mylist


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
                if type(val) == int :
                    myDict[keyA][keyB] += val
                else:
                    myDict[keyA][keyB].append(val)
            else:
                if type(val) == int :
                    myDict[keyA].update({keyB: val})
                else:
                    myDict[keyA].update({keyB: [val]})
        else:
            if type(val) == int :
                myDict.update({keyA:{keyB: val}})
            else:
                myDict.update({keyA:{keyB: [val]}})

    return myDict


def fill_double_dict_value(mydict,keylist,val):
    for g in mydict.keys():
        mydict[g]=fill_dict_value(mydict[g],keylist,val)
    return mydict


def fill_dict_value(mydict,keylist,val):
    for k in keylist:
        if k not in mydict:
            mydict[k] = val
    return mydict


def list_type(mylist,types):
    nl=[]
    for s in mylist:
        if types == "str":
            nl.append(str(s))
        elif types == "float":
            nl.append(float(s))
    return nl


# Syetem Function 
def run_command(cmd, check=True):
    """Helper function to run shell commands"""
    print(f"Running: {cmd}")
    try:
        subprocess.run(cmd, shell=True, check=check)
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {cmd}")
        print(f"   Exit code: {e.returncode}")
        if e.stderr:
            print(f"   Stderr: {e.stderr.decode() if isinstance(e.stderr, bytes) else e.stderr}")
        exit(1)


def check_path_exists(path,mode="e"):
    if not os.path.exists(path):
        if mode == "e":
            print(f"Error：Paht of software or file not exist: {path}")
            sys.exit(1) 
        elif mode == "w":
            print(f"Warming：Paht of software or file not exist: {path}")
            return False
    else:
        return True


def configure_logging():
    #logging.basicConfig(filename="EGCV.log",level=logging.INFO,format='%(asctime)s %(levelname)s %(name)s %(message)s')
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s %(name)s %(message)s',
        datefmt="%Y-%m-%d %H:%M:%S",)
    return logging.getLogger("EVCG")



## Data processing Function 

def get_sorted_chromosomes(chrom_list):
    """
    对染色体名称进行排序，支持多种格式如 lg1, lg2 或 Chr1, Chr2 等
    返回排序后的染色体列表
    """
    def extract_parts(chrom):
        # 分离字母前缀和数字部分
        chrom_str = str(chrom)
        prefix = ''.join(filter(str.isalpha, chrom_str))
        num_part = ''.join(filter(str.isdigit, chrom_str))
        
        # 处理纯数字或无数字的情况
        try:
            num = int(num_part) if num_part else 0
        except:
            num = 0
        
        return (prefix, num, chrom_str)  # 添加原始字符串作为第三元素保持稳定排序
    
    # 提取各部分信息
    chrom_info = [extract_parts(c) for c in chrom_list]
    
    # 排序：先按前缀，再按数字
    chrom_info.sort(key=lambda x: (x[0], x[1]))
    
    # 返回原始染色体名称
    return [info[2] for info in chrom_info]


def judge_significant(sign,df,column_name=None):

    cuts = []
    for cut in sign:
        if cut == "B":
            cuts.append(float(1)/len(df))
            #print(f"# The significance threshold is {cut} ")
        elif cut == "F":
            cuts.append(float(0.05)/len(df))
            #print(f"# The significance threshold is {cut}")
        elif cut == "top1":
            data = df[column_name].dropna().values
            cuts.append(np.percentile(data, 100 * (1 - 0.01)))
        elif cut == "top5":
            data = df[column_name].dropna().values
            cuts.append(np.percentile(data, 100 * (1 - 0.05)))
        else:
            try:
                cuts.append(float(cut))
            except ValueError:
                print("# -sign must be [B,F,top1,top5,float]")


    print(f"# The significance threshold is {cuts}")
    return cuts


def search_cloesd_region(gene,locis):
    # gene is a list with [chrom,start,end,strand] 
    # loci is a lsit with [[chrom,start,end,id],...]

    #print(gene)
    #print(locis[0:10])

    chrom_g, start_g, end_g, strand_g = gene[0:4]
    overlapping_regions = []
    
    for locus in locis:
        chrom_l, start_l, end_l , id_l= locus[0:4]
        
        if chrom_l != chrom_g:
            continue
        
        if (int(start_l) <= int(end_g)) and (int(end_l) >= int(start_g)):
            overlapping_regions.append(locus)

    return overlapping_regions



def merge_overlapping_intervals(df):
    """ 
    input  
        df: DataFrame, has 'Chrom', 'Start', 'End'，and sorted by Chrom and Start
    return 
        merged_df: DataFrame, in 'Chrom', 'Start', 'End' 
    """
    
    if df.empty:
        return pd.DataFrame(columns=['Chrom', 'Start', 'End'])
    
    #df['Start'] = df['Start'].astype(int)
    #df['End'] = df['End'].astype(int)
    
    merged_intervals = []
    
    current_chrom = df.iloc[0]['Chrom']
    current_start = df.iloc[0]['Start']
    current_end = df.iloc[0]['End']
    
    for i in range(1, len(df)):
        row = df.iloc[i]
        chrom = row['Chrom']
        start = row['Start']
        end = row['End']
        
        if (chrom == current_chrom and 
            start <= current_end and 
            end >= current_start):
            current_end = max(current_end, end)
        else:
            merged_intervals.append({
                'Chrom': current_chrom,
                'Start': current_start,
                'End': current_end
            })
            
            current_chrom = chrom
            current_start = start
            current_end = end
    
    merged_intervals.append({
        'Chrom': current_chrom,
        'Start': current_start,
        'End': current_end
    })
    
    merged_df = pd.DataFrame(merged_intervals)
    return merged_df


def rename_dataframe(df,chrom,pos,name,value):
    name_col= {chrom:"#CHROM",pos:"POS",value:"P",name:"ID"}
    df = df.rename(columns=name_col)
    return df


def run_parallel_thread(tasks, max_workers=4):
    """
    多线程运行函数任务（适合IO密集）
    :param tasks: [(func, args, kwargs), ...]
    :param max_workers: 线程数
    :return: results dict
    """
    results = {}

    max_workers = min(max_workers, 8)

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(func, *args, **kwargs): f"{func.__name__}{args}"
            for func, args, kwargs in tasks
        }

        for future in as_completed(futures):
            task = futures[future]
            try:
                result = future.result()
                results[task] = (0, result, None)
            except Exception as e:
                results[task] = (-1, None, str(e))

    # 👉 统一错误处理
    failed = [task for task, (code, _, _) in results.items() if code != 0]
    if failed:
        print("unfinished jobs:")
        for task in failed:
            code, _, err = results[task]
            print(f"  - {task}, exit_code={code}, error={err}")
        sys.exit(1)
    else:
        print("All jobs done!")

    return results 


def _run_cmd(cmd, capture_output):
    if capture_output:
        res = subprocess.run(cmd, shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             universal_newlines=True)  # 等价于 text=True
    else:
        res = subprocess.run(cmd, shell=True)
    return res.returncode, res.stdout if capture_output else "", res.stderr if capture_output else ""

def _run_func(func, *args, **kwargs):
    try:
        result = func(*args, **kwargs)
        return 0, result, ""
    except Exception as e:
        return 1, None, str(e)

def run_parallel2(tasks, max_workers=4, mode="cmd", capture_output=True):
    """
    并行运行命令或函数
    :param tasks: 
        - 如果 mode="cmd"，tasks 是命令字符串列表
        - 如果 mode="func"，tasks 是 (func, args, kwargs) 的列表
    :param max_workers: 最大并行数
    :param mode: "cmd" 或 "func"
    :param capture_output: 对命令是否捕获输出 (仅对 mode="cmd" 有效)
    :return: {task: (exit_code/result, stdout, stderr/exception)}
    """
    results = {}

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        if mode == "cmd":
            futures = {executor.submit(_run_cmd, cmd, capture_output): cmd for cmd in tasks}
        elif mode == "func":
            futures = {
                executor.submit(_run_func, func, *args, **kwargs): f"{func.__name__}{args}{kwargs}"
                for func, args, kwargs in tasks
            }
        else:
            raise ValueError("mode must be 'cmd' or 'func'")

        for future in as_completed(futures):
            task = futures[future]
            try:
                results[task] = future.result()
            except Exception as e:
                results[task] = (-1, None, str(e))

    failed = [task for task, (code, _, err) in results.items() if code != 0]
    if failed:
        print("unfinished jobs: ")
        for task in failed:
            code, _, err = results[task]
            print(f"  - task, exit_code={code}, error={err}")
        sys.exit(1)  
    else:
        print("All jobs done！")

    return results


if __name__ == '__main__':
    logger = configure_logging()
    main()
