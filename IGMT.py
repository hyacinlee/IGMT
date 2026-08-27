#!/usr/bin/env python3 
import argparse
import importlib
import os
import sys

#from src import Emmax
# python module need
#   numpy
#   pandas
#   seaborn
#   matplotlib
#   sklearn
#   subprocess

SUB_COMMANDS = {
    "GWAS": "GWAS", 
    "Evidence": "Evidence",
    "Phenotype": "Phenotype",
    "Visualize": "Visualize",
    "Vcftools" : "Vcftools",
    "BulkSegAna":"BulkSegAna",
    "TransAsso": "TransAsso",
    "Candidate":"Candidate",
    "GeneEffect":"GeneEffect",
    "Haplotype":"Haplotype"
}

def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description='''
    Main Prgram of Multi-omics genetic mapping tools , Plase choose Function from:

        #--------------------------------  Main module   --------------------------------#

        Vcftools     :   The module for processing vcf file , including format conversion, adding or simplifying information, renaming chromosomes or sample names, etc.

        Phenotype    :   The module for processing phenotypic data, including clustering, dimensionality reduction, etc. 

        GWAS         :   The module for conducting association analysis using Emmax.

        Evidence     :   The module for evidence genes use multi-omics results.

        Visualize    :   The module for output Manhantun or other figures. 

        BulkSegAna   :   The module for perform Bulked Segregant Analysis use vcf files.


        #------------------------------  Additional module  ------------------------------#

        GeneEffect  :    The module for Post-GWAS to identify gene-level effect.

        Candidate   :    The module for select candidate genes using muti-methods.

        Haplotype    :   The module for perform Haplotype analysis based on single gene.

        TransAsso    :   The module for perform Association analysis between gene expression matrix and phenotype in population.
        
    Date      : 2026-08-25
    Version   : 1.0.0
        ''',
        add_help=False,
        usage="%(prog)s < Function > [args]"
    )
    parser.add_argument(
        "Function",
        choices=SUB_COMMANDS.keys(),
        help="Avaliable Functions",
    )
   
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    

    # add src to pathon libs 
    script_dir = os.path.dirname(os.path.abspath(__file__))
    src_dir = os.path.join(script_dir, 'src')
    if src_dir not in sys.path:
        sys.path.append(src_dir)

 
    if "-h" in sys.argv or "--help" in sys.argv:
        if len(sys.argv) == 2:
            parser.print_help()
            sys.exit(0)
        else:
            cmd = sys.argv[1]
            main_func = load_module(cmd)
            main_func(["-h"])  
            sys.exit(0)

    # 正常解析命令名
    args, remaining_args = parser.parse_known_args()
    
    # 加载并执行子命令
    main_func = load_module(args.Function)
    main_func(remaining_args) # 传递剩余参数

def load_module(module_name):
    """动态加载子命令模块并返回其 main 函数"""
    try:
        module = importlib.import_module(SUB_COMMANDS[module_name])
        return module.main
    except ImportError as e1:
        print(f"ImportError:\n{e1}")
        sys.exit(1)
    except KeyError as e2:
        print(f"KeyError:\n{e2}")
        sys.exit(1)


if __name__ == "__main__":
    main()
