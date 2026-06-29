# IGMT

IGMT (Integrative Genetic Mapping Toolkit), an integrative framework for large-scale phenotypic GWAS that unifies phenotype imputation, trait clustering, association testing, and gene prioritization. IGMT leverages machine-learning models for imputation and trait complexity reduction, and integrates multi-omics and genetic evidence to improve candidate gene discovery, ultimately enabling efficient and interpretable extraction of biological insights from large-scale phenotypic data.

## Install
 This toolkit is entirely developed based on python. 
 Please ensure that your python >= 3.6 and pip3 in your env:


Then :
```

cd ${your_path}
git clone https://github.com/hyacinlee/IGMT

cd ./IGMT
pip3 install -r requirements.txt

export PATH="${your_path}/IGMT/:$PATH"
export PYTHONPATH="${your_path}/IGMT/src:$PYTHONPATH"

IGMT.py -h

```


## Module Introduction
* **Vcftools**          :   The module for processing vcf file , including format conversion, adding or simplifying information, renaming chromosomes or sample names, etc.
* **Phenotype**         :   The module for processing phenotypic data, including clustering, dimensionality reduction, etc. 
* **GWAS**              :   The module for conducting association analysis using Emmax.
* **Evidence**          :   The module for evidence genes use multi-omics results.
* **Visualize**         :   The module for output Manhantun or other figures.
* **GeneEffect**        :   The module for post-GWAS gene-level effect rank.
* **Candidate**         :   The module for select candidate genes using muti-methods.
* **BulkSegAna**        :   The module for perform Bulked Segregant Analysis use vcf files.
* **TransAsso**         :   The module for perform Association analysis between gene expression matrix and phenotype in population.


## Vcftools
```
# rename vcf chromosome names 
IGMT.py Vcftools -v w220RNA.PASS.filter.strict.imputed.vcf --rename_chr --chr TGY.rename -o TGY.final.vcf
```


##  Phenotype 

#### Basic matrix operations
```
Extract subset (by row or column):
------------------------------------------------------------
IGMT.py Phenotype -p example.raw.phe --subLines Candidate.triats.ids --out Candidate.triats
IGMT.py Phenotype -p example.raw.phe --subColumns Candidate.sample.ids --out Candidate.sample


Run statistics by row:
------------------------------------------------------------
IGMT.py Phenotype -p example.raw.phe --stat -o example


Combine "combine_files" into the base file (-p):
------------------------------------------------------------
IGMT.py Phenotype --combine Sample --combine_method Base -p example.raw.phe --combine_files new.traits1.txt new.traits2.txt  --out Combine


Column-wise group mean calculation (typically used for handling multiple replicates):
------------------------------------------------------------------------------------------------------------------------
IGMT.py Phenotype -p gene.expression.phe --mean group.txt --out  gene.expression 
head group.txt 
colum_name1 group1
colum_name2 group1
colum_name3 group2
colum_name4 group2
...

# t-test
IGMT.py Phenotype -p gene.expression.grouped.tsv --group_ttest t-group.txt --out gene.expression
head t-group.txt
group1  Hight
group2  Hight
group3  Low
group4  Low
...

```

#### Impution
```
# Impute missing values using the input matrix itself
IGMT.py Phenotype -p Area_matrix_final_filter.tsv --impute iterative_rf --top_k 10 --corr_threshold 0.7 -o Area_matrix_final_filter_impute1.tsv

# Impute missing values using an external database matrix
IGMT.py Phenotype -p Area_matrix_final_filter_impute1.tsv  --impute iterative_rf --top_k 10 --corr_threshold 0.7 --database ../Date_base.met3971.rename.tsv -o Area_matrix_final_filter_impute2.tsv

```


#### Clustering
```
Global clustering with 500 iterations
------------------------------------------------------------
IGMT.py Phenotype -p ../example.raw.phe -g -c 0.7 -z 0.8 -o meta.know.C70.Z80 -b 500


Use the specified trait for local clustering (search all traits for those similar to `-l`)
------------------------------------------------------------
IGMT.py Phenotype -p ../example.raw.phe -l local.txt -c 0.7 -z 0.8 -o meta.know.C70.Z80.local -b 20
head local.txt 
YL_cmp1651
YL_cmp1307
YL_cmp624
YL_cmp1304
YL_cmp622
...

```


---
## GWAS 
```
# Basic GWAS
IGMT.py Emmax -v ./snp.dep5.mis90.maf5.simply.vcf -c PCA.cov -p phe.txt

# If there are too many SNPs (> 2000,000), use parallel processing to improve efficiency
IGMT.py Emmax -v ./snp.dep5.mis90.maf5.simply.vcf -c PCA.cov -p phe.txt --split chrom.lst

# After completion, adjust the threshold and re-filter the results
rm -rf GWAS.outputs.done  GWAS.manhantun.done
IGMT.py Emmax -v ./snp.dep5.mis90.maf5.simply.vcf -c PCA.cov -p phe.txt --split chrom.lst -s 0.000001 

# Sample File (The order is not important)
head PCA.cov 
Sample  PC1 PC2 PC3
10  0.067394    0.0363772   0.023049
104 0.0326593   0.0115734   -0.0513225
105 0.0652511   0.061142    0.0384592
108 0.0467105   0.0194756   -0.017885
110 0.0614426   0.0341291   0.0129854

head chrom.lst
lg1
lg2
lg3
lg4
lg5
```



---
## Evidence
cmd:
```
IGMT.py Evidence -i ../03.Gwas_singleTriat/Traits.Pos2.emmax.ps.significant.tsv -p ../03.Gwas_singleTriat/Traits.Pos2.phe -o Traits.Pos2.200k -a ../../../../../02.Genotype/02.Genotype_consGenome/03.Annovar/snp.dep5.mis90.maf5.simply -f /work/minghui/Tea_reseq/3.Reseq_FDDB/00.data/refgenome/ref_one/Function.swissport.txt -b /work/minghui/Tea_reseq/3.Reseq_FDDB/00.data/refgenome/ref_one/FD.mRNA.bed -v ../../../4.Area.Grouped.cons/data/snp.dep5.mis90.maf5.simply.vcf.gz --regions 1 200 --regions_score 2 1  -e ../../../../00.data/Evidence/filter.eqtl.cis.txt -w ../../../../00.data/Evidence/weight.txt
```

#### Weight file
```
#Name   Type    Weight  File 
eQTL    eqtl    5       ./add/187_filtered_exp_matrix_qqnorm_file_sort_snp_matrixEQTLoutput_signE5.txt
PopFpkm express 5       ./add/DASZ.220.exp.sorted.FPKM.head.txt
SVs     loci    3       ./add/More2.cluster.svs.bed
trans1  gene    2       ./add/low_vs_high.DEGlist.txt
YL_cmp781       gene    1       ./add/Traits.YL_cmp781.sign.BasicGene.xls
YL_cmp828       gene    1       ./add/Traits.YL_cmp828.sign.BasicGene.xls
```

---
## Visualize

#### plot manhantun
```
IGMT.py Visualize -m manha --site_value ../03.Gwas_singleTriat/Traits.Pos2.emmax.ps.tsv --hightlight candidate.genes.heightlight --out Traits.Pos2.heightlight


cat heightlight file:
lg12_16125300   PPOX2   red
lg12_35022415   Cyp450  red
lg13_193074629  LAP2    red
lg3_83862148    AAP3    red
```


#### plot boxplot
```
IGMT.py Visualize -m stat --df Test.group.combine.tsv --box -o Test.group.combine --xData Group --yData Pos2 --ttest CSA1,CSS1 CSA1,CSS2 CSA2,CSS2 CSA2,CSS1 --figsize 8 8 --group_color Group.color

cat Group.color 
CSS1    #3274a1
CSS2    #e1812c
CSA1    #FB5607
CSA2    #8338EC

```


#### plot Cor
```
IGMT.py Visualize -m stat --df Test.group.combine.tsv --cor -o Test.group.combine --figsize 8 8 --exc Group --plot Row
```



---
### Candidate
```

cat Coloc.txt
Pos2    ../../03.Gwas_singleTriat/Traits.Pos2.500k.sign.WeightGenes.xls
Pos5    ../../03.Gwas_singleTriat/Traits.Pos5.500k.sign.WeightGenes.xls
Pos15   ../../03.Gwas_singleTriat/Traits.Pos15.500k.sign.WeightGenes.xls
Pos35   ../../03.Gwas_singleTriat/Traits.Pos35.500k.sign.WeightGenes.xls
Pos53   ../../03.Gwas_singleTriat/Traits.Pos53.500k.sign.WeightGenes.xls

IGMT.py Candidate  -i Coloc.txt --co_local  -o Trait27.express --fit_exp "G128!=None"

```



---
## BSA
```
# BSA in single SNP site:
IGMT.py BulkSegAna -v snp.filter.vcf -s1 Low -s2 High -t ed -s top1 top5 

# BSA in a 50kb-step and 500kb-windows sliding windows:
IGMT.py BulkSegAna -v snp.filter.vcf -s1 Low -s2 High -t ed -s top1  -step 50  -wind 500 -minSite 5  -minDepth 10


```



---
## Trans Associate
```
# Although the covariance file (-c) is not mandatory, we strongly recommend using it.

# Only Run two traits:
IGMT.py TransAsso -e DASZ.220.exp.filter1.sub_set.tsv -p ../GWAS.phe  -t Group1 Group2   -c ./PCA.matrix

# Run all traits in .phe file with parallel threads:
IGMT.py TransAsso -e DASZ.220.exp.filter1.sub_set.tsv -p ../GWAS.phe --threads 24  -c ./PCA.matrix

# Fig out candidate genes (Only the gene ID in the first column is needed)
IGMT.py TransAsso -e DASZ.220.exp.filter1.sub_set.tsv -p ../GWAS.phe  -t Group3 -c ./PCA.matrix --figGenes TransAsso.Group3.signal.txt
```


---
### License

This project is licensed under the MIT License.
See the [LICENSE](LICENSE) file for details.

---
### Citation

