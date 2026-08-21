## Introduction

IGMT provides a set of modules covering phenotype preprocessing, phenotype imputation and clustering, genome-wide association analysis, candidate gene prioritization, gene-level effect estimation, haplotype analysis, BSA, transcriptome–phenotype association, VCF processing, and visualization.

The major modules are:

- **Vcftools**: basic processing of VCF files, including filtering, simplification, and renaming of chromosomes or samples.
    
- **Phenotype**: phenotype matrix processing, missing-value imputation, trait correlation analysis, and trait clustering.
    
- **GWAS**: multi-trait genome-wide association analysis based on EMMAX.
    
- **Evidence**: candidate gene prioritization by integrating GWAS signals, variant annotation, genomic position, eQTLs, and additional multi-omics evidence.
    
- **GeneEffect**: gene-level effect estimation from SNP genotypes using ElasticNet or PCA-based aggregation.
    
- **Candidate**: integration and co-localization analysis of candidate genes obtained from multiple traits or analyses.
    
- **Haplotype**: haplotype analysis for specified candidate genes and comparison of haplotypes with phenotypic variation.
    
- **BulkSegAna**: bulked segregant analysis using Euclidean distance or SNP-index approaches.
    
- **TransAsso**: population-level association analysis between gene expression and phenotypes.
    
- **Visualize**: visualization of GWAS results, phenotype distributions, correlations, PCA, and genotype–phenotype relationships.
    

For the complete list of command-line options for each module, use:

```bash
IGMT.py <Module> -h
```

## Install

This toolkit is entirely developed based on python. Please ensure that your python >= 3.6 and pip3 in your env:

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
---

# 1. Vcftools

The `Vcftools` module provides lightweight VCF preprocessing functions required by downstream IGMT analyses. It is intended for common format-processing operations rather than as a replacement for specialized VCF-processing software.

### Example

```bash
# Rename chromosomes
IGMT.py Vcftools \
    -v input.vcf \
    --rename_chr \
    --chr chromosome.rename.txt \
    -o output.rename.vcf

# add snp id infomation if
IGMT.py Vcftools \
    -v input.rename.vcf \
    --add_ids \
    -o output.add.vcf
```

### Main parameters

|Parameter|Description|
|---|---|
|`-v, --vcf`|Input VCF file|
|`-c, --chr`|Chromosome list or chromosome-renaming file|
|`-k, --keep`|Sample list or sample-renaming file|
|`-t, --tags`|FORMAT fields to retain, such as `GT AD DP`|
|`-o, --out`|Outfile name|
|`-a,--add_ids`|Generate SNP IDs in `chromosome_position` format|
|`--simplify`|Remove unnecessary VCF information and retain genotype information|
|`--rename_chr`|Rename chromosomes according to `--chr`|
|`--rename_sample`|Rename samples according to `--keep`|
|`--diploid`|Convert haploid genotypes to pseudo-diploid genotypes|

### Parameter selection

These parameters mainly control file formatting and therefore do not affect statistical inference. Users should retain the FORMAT fields required by subsequent analyses and ensure that chromosome and sample identifiers are consistent among genotype, phenotype, covariate, and annotation files.

### Result interpretation

The output is a processed VCF file. Before downstream analysis, users should confirm that chromosome names, SNP IDs, sample IDs, and genotype FORMAT fields are consistent with the corresponding phenotype and annotation files.

---

# 2. Phenotype

The `Phenotype` module performs phenotype matrix preprocessing, missing-value imputation, correlation analysis, and trait clustering. It is designed particularly for large phenotypic datasets containing many correlated traits.

## 2.1 Basic matrix operations

```bash
# Extract traits
IGMT.py Phenotype \
    -p example.raw.phe \
    --subLines Candidate.traits.ids \
    -o Candidate.traits

# Extract samples
IGMT.py Phenotype \
    -p example.raw.phe \
    --subColumns Candidate.sample.ids \
    -o Candidate.samples

# Calculate row-wise statistics
IGMT.py Phenotype \
    -p example.raw.phe \
    --stat \
    -o example
```

The module also supports matrix merging, transposition, replicate averaging, correlation calculation, and group-based t-tests.

---

## 2.2 Missing-value imputation

IGMT implements two iterative imputation strategies:

- `iterative_rf`: iterative random-forest regression
    
- `iterative_bayes`: iterative Bayesian regression
    

For each phenotype with missing observations, other traits correlated with the target phenotype are selected as predictors. Predictor selection is controlled by the correlation threshold and the maximum number of predictors.

### Example

```bash
IGMT.py Phenotype \
    -p phenotype.tsv \
    --impute iterative_rf \
    --top_k 10 \
    --corr_threshold 0.7 \
    --n_iter 5 \
    -o phenotype.imputed.tsv
```

An external phenotype or multi-dimensional reference matrix can also be supplied:

```bash
IGMT.py Phenotype \
    -p phenotype.tsv \
    --impute iterative_rf \
    --database reference.matrix.tsv \
    --top_k 10 \
    --corr_threshold 0.7 \
    -o phenotype.imputed.tsv
```

### Important parameters

|Parameter|Default|Description|
|---|--:|---|
|`--impute`|—|`iterative_rf` or `iterative_bayes`|
|`--corr_threshold`|0.6|Minimum correlation required for a trait to be considered as a predictor|
|`--top_k`|10|Maximum number of the most strongly correlated predictor traits|
|`--n_iter`|5|Maximum number of iterative imputation rounds|
|`--threads`|24|Number of parallel threads|
|`--database`|None|Optional external matrix used to provide predictor traits|

### Parameter selection and expected effects

`--corr_threshold` controls how strongly related a predictor must be to the target phenotype. A higher threshold selects fewer but more strongly correlated traits, whereas a lower threshold increases the number of available predictors but may introduce weakly informative variables.

`--top_k` limits the maximum number of correlated traits used for prediction. Increasing this value provides more predictors when many related traits are available, but does not necessarily improve imputation when additional traits are weakly related.

`--n_iter` controls the maximum number of iterative updating rounds. The default value is suitable for routine analyses; increasing it may be useful when missingness is extensive but also increases computation time.

Because correlation structures differ substantially among datasets and traits, the default parameters should be regarded as practical starting values rather than universally optimal settings. For important datasets, users are encouraged to evaluate several combinations of `--corr_threshold` and `--top_k`, for example by artificially masking a subset of observed values and comparing imputation accuracy.

### Result interpretation

The main output is the completed phenotype matrix. Imputed values should not automatically be considered equivalent to directly measured observations. For downstream analyses, users should examine the proportion of missing values, the number and strength of available predictor traits, and, where possible, evaluate imputation accuracy using masked observed values.

Traits with few or no sufficiently correlated predictors may be less reliably imputed than traits embedded in strongly correlated phenotype groups.

---

## 2.3 Trait clustering

IGMT implements iterative correlation–principal component clustering to group highly correlated phenotypes and reduce redundant traits in large-scale phenotype datasets.

### Global clustering

```bash
IGMT.py Phenotype \
    -p example.raw.phe \
    -g \
    -c 0.7 \
    -z 0.8 \
    -b 500 \
    -o meta.C70.Z80
```

### Local clustering

```bash
IGMT.py Phenotype \
    -p example.raw.phe \
    -l local.txt \
    -c 0.7 \
    -z 0.8 \
    -b 500 \
    -o meta.C70.Z80.local
```

### Important parameters

|Parameter|Default|Description|
|---|--:|---|
|`-g, --globals`|—|Perform clustering for all traits|
|`-l, --locals`|—|Perform clustering around specified traits|
|`-c, --corr_threshold`|0.6|Minimum pairwise correlation used when constructing a candidate group|
|`-z, --pc_threshold`|0.8|Minimum PC1 contribution required for retaining traits in a cluster|
|`-b, --bootstrap`|500|Number of repeated clustering iterations|
|`-f, --fig`|False|Output correlation heatmaps for resulting groups|

### Parameter selection and expected effects

The correlation threshold (`-c`) controls the initial similarity required among traits. Increasing this threshold generally produces smaller and more homogeneous groups, whereas decreasing it allows more moderately correlated traits to enter the same candidate group.

The PC1 threshold (`-z`) controls the amount of common variation that must be captured by the first principal component. A higher value imposes a stricter requirement for within-group consistency.

Because the procedure begins from iteratively selected base traits, `-b` repeats the clustering process and reduces uncertainty associated with a particular starting order. The default of 500 iterations is recommended for final analyses; fewer iterations can be used during exploratory testing.

### Result interpretation

Each resulting group represents a set of phenotypes sharing substantial correlation and a common major variation component. A group containing many traits indicates a strong redundant or coordinated phenotype pattern, whereas isolated traits generally represent relatively independent phenotypic variation.

The derived principal-component phenotype can subsequently be used as a representative phenotype in GWAS, reducing the number of highly redundant association tests.

---

# 3. GWAS

The `GWAS` module performs genome-wide association analysis for one or multiple phenotypes using EMMAX.

### Basic analysis

```bash
IGMT.py GWAS \
    -v snp.vcf \
    -p phenotype.phe \
    -c PCA.cov
```

For large VCF files, chromosomes can be processed separately:

```bash
IGMT.py GWAS \
    -v snp.vcf \
    -p phenotype.phe \
    -c PCA.cov \
    --split chrom.lst \
    --threads 24
```

### Important parameters

|Parameter|Default|Description|
|---|--:|---|
|`-v, --vcf`|—|Input genotype VCF|
|`-p, --phe`|—|Phenotype file|
|`-c, --cov`|None|Covariates such as principal components or population structure|
|`-m, --maf`|0.05|Minimum minor allele frequency|
|`-s, --sign`|`B`|Significance threshold|
|`-t, --threads`|12|Number of parallel threads|
|`--split`|None|Chromosome list for chromosome-wise parallel analysis|

For `--sign`:

- `B`: `1 / number_of_SNPs`
    
- `F`: `0.05 / number_of_SNPs`
    
- a numeric value: user-defined P-value threshold
    

### Parameter selection and expected effects

The MAF threshold removes rare variants with insufficient genotype counts for reliable association testing. Lowering the threshold retains more rare variants but may reduce statistical stability, particularly in populations with modest sample sizes.

Population structure or other relevant covariates should normally be included using `-c` to reduce confounding.

`--split` changes computational efficiency rather than the statistical model and is recommended for very large genotype datasets.

### Main outputs and interpretation

For each phenotype, IGMT generates GWAS result tables and Manhattan plots. The result table contains SNP positions and association statistics.

Lower P-values indicate stronger evidence for phenotype–genotype association. SNPs exceeding the selected significance threshold are written to significant-site output files and can be passed directly to the `Evidence`, `GeneEffect`, or `Haplotype` modules.

Genome-wide significance should be interpreted together with local LD structure, allele frequency, population structure, and biological annotation rather than solely from the most significant SNP.

---

# 4. Evidence

The `Evidence` module prioritizes candidate genes by integrating significant GWAS SNPs with variant functional annotation, physical distance to genes, eQTL information, and optional additional biological evidence.

### Example

```bash
IGMT.py Evidence \
    -i Traits.Pos2.emmax.ps.significant.tsv \
    -v genotype.vcf.gz \
    -a snp.annovar \
    -b gene.bed \
    -p Traits.Pos2.phe \
    -e cis_eQTL.txt \
    -w weight.txt \
    --regions 1 200 \
    --regions_score 2 1 \
    -o Traits.Pos2
```

### Core evidence parameters

|Parameter|Default|Description|
|---|--:|---|
|`--snp_scores`|`10 5 3 3`|Weights for nonsynonymous, synonymous, intronic, and splicing SNP categories|
|`--regions`|`1`|Distances from associated SNPs to nearby genes, in kb|
|`--regions_score`|`1`|Weight assigned to each distance category|
|`--eQTL_score`|3|Weight assigned to SNP–gene eQTL evidence|
|`--base_scores`|0.8|Relative contribution of genetic evidence to the final integrated score|
|`--min_sites`|5|Minimum number of candidate SNPs required to continue analysis|
|`-w, --weight`|None|Optional additional evidence and corresponding weights|

The number of values supplied to `--regions_score` must match the number supplied to `--regions`.

### Interpretation of the default weights

The default SNP-category weights represent the expected relative functional importance of different forms of genetic evidence rather than statistically estimated effect sizes. Therefore, nonsynonymous variants receive a larger default contribution than synonymous, intronic, or nearby variants.

These values can be modified according to the biological context and the type of available data.

### Additional evidence

Additional evidence can be supplied using a weight file:

```text
#Name       Type       Weight    File
Expression  express    5         expression.tsv
SV          loci       3         structural_variants.bed
DEG         gene       2         DEG.list
```

Supported evidence types include:

- `express`: gene-expression evidence
    
- `loci`: genomic interval evidence
    
- `gene`: candidate-gene list evidence
    

The relative contributions of individual additional evidence types can be modified using the `Weight` column.

### Combining genetic and additional evidence

When additional evidence is supplied, IGMT calculates:

```text
Total_score =
    base_scores × Base_score
    + (1 - base_scores) × Evi_score
```

The default `--base_scores 0.8` therefore gives greater weight to the original genetic evidence while allowing complementary multi-omics information to modify candidate-gene ranking.

Increasing `--base_scores` makes the ranking more dependent on GWAS and variant evidence; decreasing it increases the influence of additional biological evidence.

### Main outputs

Typical outputs include:

```text
*.sign.BasicSite.xls
*.sign.BasicGene.xls
*.sign.WeightGenes.xls
```

`BasicSite` summarizes associated SNPs, their functional classes, mapped genes, distance information, and GWAS statistics.

`BasicGene` summarizes SNP-based evidence at the gene level and reports a `Base_score`.

When additional evidence is provided, `WeightGenes` reports the integrated `Total_score`, together with genetic and additional evidence scores.

### Result interpretation

Genes with larger `Base_score` have stronger support from the underlying GWAS-associated variants and their functional categories.

Genes with larger `Total_score` receive stronger combined support after incorporating additional evidence. Candidate genes should normally be interpreted by considering both the final ranking and the specific evidence contributing to that ranking rather than treating the score as a probability of causality.

---

# 5. GeneEffect

The `GeneEffect` module estimates the phenotypic contribution of SNPs assigned to each gene. SNPs are first filtered according to LD and are then summarized using either ElasticNet regression or PCA before estimating gene-level phenotype variance.

### Example

```bash
IGMT.py GeneEffect \
    -s SNP_gene.annotation.tsv \
    -g genotype.tsv \
    -p phenotype.phe \
    -c PCA.cov \
    --method ElasticNet \
    --ld-r2 0.9 \
    --min-snps 1 \
    -o Trait1
```

### Important parameters

|Parameter|Default|Description|
|---|--:|---|
|`--method`|ElasticNet|SNP aggregation method: `ElasticNet` or `PCA`|
|`--ld-r2`|0.9|LD threshold used for SNP pruning within genes|
|`--min-snps`|1|Minimum number of SNPs required per gene after pruning|
|`--l1-ratio`|0.5|ElasticNet mixing parameter|
|`-c, --cov`|None|Optional covariates|

### Parameter selection and expected effects

A lower `--ld-r2` applies more stringent LD pruning and reduces redundancy among SNPs, whereas a higher threshold retains more correlated SNPs.

`ElasticNet` is suitable when multiple SNPs may jointly contribute to the phenotype and provides regularized SNP weighting. `PCA` provides a simpler representation of the major genotype variation within each gene.

The default `--l1-ratio 0.5` balances L1 and L2 regularization. Users generally do not need to modify this value unless they specifically wish to change the sparsity of SNP effects.

### Result interpretation

The primary output is:

```text
*.gene_effect.ElasticNet.tsv
```

or the corresponding PCA output.

Genes with larger gene-level explained-effect statistics contribute more strongly to phenotypic variation in the analyzed population. These estimates provide a complementary ranking to SNP-level GWAS significance and the evidence-based ranking produced by the `Evidence` module.

A high gene-level effect does not by itself establish causality and should be interpreted together with GWAS significance and biological evidence.

---

# 6. Candidate

The `Candidate` module integrates candidate-gene results across multiple traits and identifies genes repeatedly supported by related phenotypes.

### Input file

```text
Pos2    Traits.Pos2.sign.WeightGenes.xls
Pos5    Traits.Pos5.sign.WeightGenes.xls
Pos15   Traits.Pos15.sign.WeightGenes.xls
```

### Example

```bash
IGMT.py Candidate \
    -i Coloc.txt \
    --co_local \
    --fit_exp "G128!=None" \
    -o TraitGroup
```

### Main parameters

|Parameter|Description|
|---|---|
|`--co_local`|Integrate candidate genes across multiple traits|
|`--fit_exp`|Filter input candidate tables using a specified column expression|
|`-o`|Output prefix|

### Result interpretation

The module reports both the merged candidate-gene records and a gene-level summary across traits.

Genes supported by multiple related traits may represent shared genetic regulators or pleiotropic candidate genes. However, repeated detection among correlated phenotypes should not automatically be interpreted as independent evidence; the result is intended primarily for prioritization.

---

# 7. Haplotype

The `Haplotype` module performs candidate-gene-level haplotype analysis using GWAS SNPs and genotype data and relates the resulting haplotypes to phenotypic values.

### Example

```bash
IGMT.py Haplotype \
    -i Traits.Pos2.emmax.ps.tsv \
    -g CandidateGene \
    -p Traits.Pos2.phe \
    -v genotype.vcf.gz \
    -a snp.annovar \
    -b genes.bed \
    -r NoSyn CDS \
    -o CandidateGene
```

### Important parameters

|Parameter|Default|Description|
|---|--:|---|
|`-g, --gene`|—|Candidate gene|
|`-r, --region`|All|SNP classes included: `All`, `CDS`, `NoSyn`, `intronic`, or `Promoter`|
|`-m, --method`|Hom|Haplotype processing method|
|`-f, --flank`|2000|Promoter flanking distance in bp|

### Parameter selection

`--region` determines which variants are used to define haplotypes. Restricting the analysis to coding or nonsynonymous variants provides a more functionally focused analysis, while `All` retains the maximum available genetic information.

### Main outputs and interpretation

Typical outputs include:

```text
*.Site_info.txt
*.hap_info.txt
*.sample_hap.txt
```

`Site_info` describes the variants used to construct haplotypes, `hap_info` summarizes haplotype groups and their phenotype information, and `sample_hap` reports the haplotype assigned to each sample.

Phenotypic differences among haplotypes provide supporting evidence that local genetic variation around a candidate gene is associated with the trait.

---

# 8. BulkSegAna

The `BulkSegAna` module performs bulked segregant analysis using either Euclidean distance (`ed`) or SNP-index (`index`).

### Single-site analysis

```bash
IGMT.py BulkSegAna \
    -v snp.filter.vcf \
    -s1 Low \
    -s2 High \
    -t ed \
    -s top1
```

### Sliding-window analysis

```bash
IGMT.py BulkSegAna \
    -v snp.filter.vcf \
    -s1 Low \
    -s2 High \
    -t ed \
    -step 50 \
    -wind 500 \
    -minSite 5 \
    -minDepth 10
```

### Important parameters

|Parameter|Default|Description|
|---|--:|---|
|`-t`|ed|`ed` or `index`|
|`-step`|None|Sliding-window step size in kb|
|`-wind`|None|Sliding-window size in kb|
|`-minDepth`|10|Minimum sequencing depth in each bulk|
|`-minSite`|10|Minimum number of SNPs per window|
|`-s`|top1|Candidate-region threshold|

Smaller windows provide higher spatial resolution but contain fewer SNPs and may therefore be noisier. Larger windows provide smoother signals but reduce localization resolution.

Candidate regions are represented by genomic intervals showing unusually large differences between the two bulks.

---

# 9. TransAsso

The `TransAsso` module tests associations between population-level gene-expression variation and phenotypic variation using linear regression.

### Example

```bash
IGMT.py TransAsso \
    -e expression.tsv \
    -p phenotype.phe \
    -c PCA.matrix \
    --threads 24
```

Specific traits can be selected using:

```bash
IGMT.py TransAsso \
    -e expression.tsv \
    -p phenotype.phe \
    -t Trait1 Trait2 \
    -c PCA.matrix
```

### Important parameters

|Parameter|Description|
|---|---|
|`-e, --express`|Gene-expression matrix|
|`-p, --phe`|Phenotype matrix|
|`-c, --cov`|Covariates|
|`-t, --trait`|Traits to analyze; all traits are analyzed if omitted|
|`-l, --loci`|Optional gene-position file|
|`--threads`|Parallel threads|
|`--figGenes`|Genes for which expression–phenotype plots should be generated|

### Parameter selection

Although covariates are optional, population-structure or other appropriate covariates are strongly recommended when expression and phenotype data are measured in a structured population.

### Result interpretation

Genes with small association P-values show expression levels that covary with the phenotype after accounting for the supplied covariates.

Such associations provide transcriptomic support for candidate genes but do not distinguish whether expression variation is causal, consequential, or correlated with another underlying factor.

---

# 10. Visualize

The `Visualize` module provides common visualization functions used throughout IGMT.

Three major modes are available:

```text
manha    Manhattan plots
pvg      genotype–phenotype visualization
stat     general statistical visualization
```

### Manhattan plot

```bash
IGMT.py Visualize \
    -m manha \
    --site_value Traits.Pos2.emmax.ps.tsv \
    --hightlight candidate.genes.highlight \
    -o Traits.Pos2
```

Important options include `--sign`, `--sub`, and `--local`, which control significance thresholds and genomic regions displayed.

### Boxplot

```bash
IGMT.py Visualize \
    -m stat \
    --df Test.group.combine.tsv \
    --box \
    --xData Group \
    --yData Pos2 \
    --ttest CSA1,CSS1 CSA1,CSS2 \
    --figsize 8 8 \
    -o Test.group.combine
```

### Correlation heatmap

```bash
IGMT.py Visualize \
    -m stat \
    --df Test.group.combine.tsv \
    --cor \
    --plot Row \
    --exc Group \
    -o Test.group.combine
```

### Main visualization parameters

|Parameter|Default|Description|
|---|--:|---|
|`--figsize`|16 5|Figure dimensions|
|`--dpi`|300|PNG resolution|
|`--ftype`|png|`png` or `pdf`|
|`--sign`|B|Manhattan-plot significance threshold|
|`--local`|None|Plot a specified genomic interval|
|`--hightlight`|None|Highlight specified SNPs|
|`--cor`|False|Correlation heatmap|
|`--pca`|False|PCA visualization|
|`--box`|False|Boxplot|
|`--lm_dot`|False|Linear-regression scatter plot|

Visualization parameters affect presentation rather than the statistical analyses themselves.

---

# Recommended workflow

A typical large-scale phenotypic GWAS workflow in IGMT is:

```text
Phenotype preprocessing
        ↓
Missing-value imputation
        ↓
Trait correlation / clustering
        ↓
GWAS
        ↓
Evidence-based candidate gene prioritization
        ↓
GeneEffect / Candidate
        ↓
Haplotype or TransAsso validation
        ↓
Visualization
```

Not every module is required for every dataset. Individual modules can be used independently according to the available data and the biological question.

## Parameter considerations

The default parameters distributed with IGMT are intended to provide practical starting points for routine analyses. They should not be considered universally optimal for all datasets.

In particular:

- phenotype imputation parameters should be selected according to trait correlation structure and missing-data patterns;
    
- trait-clustering thresholds determine the balance between cluster size and within-cluster consistency;
    
- GWAS thresholds and covariates should reflect marker density and population structure;
    
- candidate-gene evidence weights represent user-adjustable biological priorities rather than statistically estimated effect sizes;
    
- additional omics evidence should only be included when the corresponding data are available and biologically relevant.
    

For important analyses, users are encouraged to evaluate parameter sensitivity and retain the parameter settings used for each analysis to facilitate reproducibility.
