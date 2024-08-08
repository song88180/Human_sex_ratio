# In search of the genetic variants of human sex ratio at birth: Was Fisher wrong about sex ratio evolution?

Here are analysis scripts used in (Song & Zhang, 2024, under review) titled 
"In search of the genetic variants of human sex ratio at birth: Was Fisher wrong about sex ratio evolution?"

### Required software and packages:

The versions of software/packages I used are shown in the parentheses. You 
can of course use different versions that suit your needs

* Python (v3.10.13)
  * NumPy (v1.26.2)
  * pandas (v2.1.1)
  * SciPy (v1.10.1)
  * Matplotlib (v3.7.1)
  * seaborn (v0.13.2)
  * scikit-learn (v1.2.2)
  * statsmodels (v0.13.5)
  * pickle (v4.0)
* R (v4.3.1)
  * sumFREGAT (v1.2.5)
* SLiM (v3.7)
* Regenie (v3.4)

### 01\_GWAS\_OSR\_cov\_logistic.py
The script performs GWAS on siblings' sex to find potential SNPs underlying human offspring sex ratio. The results are shown in **Figure 2a** of the paper.

### 02.1\_UKB\_power\_analysis.ipynb
The script performs statistical power analysis for the GWAS on siblings' sex by simulation. The results are shown in **Figure 1c** and **Figure S1**of the paper.

### 02.2\_Zietsch\_power\_analysis.ipynb
The script performs statistical power analysis for the method in Zietsch *et al.*'s study. The results are shown in **Figure 1a** of the paper.

### 02.3\_Boraska\_power\_analysis.ipynb
The script performs statistical power analysis for the method in Boraska *et al.*'s study. The results are shown in **Figure 1a** of the paper.

### 03\_Odds\_ratio\_correction.ipynb
The script converts the odds ratio of M/F in siblings to odds ratio of M/F in offsprings by simulation. The results are shown in **Figure S2** of the paper.

### 04\_Family\_based\_heritability.ipynb
The script calculates the family-based heritability of the offspring sex ratio by using relative pairs data from the UK Biobank. The results are shown in **Figure 3b** of the paper.

### 05.1\_Gene\_test\_sumFREGAT.ipynb
The script utilizes the R package "sumFREGAT" to conduct gene-based test with generated GWAS summary statistics. The results are shown in **Figure S3a** of the paper.

### 05.2\_burden\_test\_write\_mask\_bed.sh
The script generates a plink bed file that contained burden scores of all genes for all UKB individuals using the “--write-mask” option in Regenie.

### 05.3\_burden\_test\_GWAS\_cov\_logistic.py
The script performs genome-wide gene burden test on siblings' sex to find potential genes influencing human offspring sex ratio. The results are shown in **Figure S3b** of the paper.

### 06.1_Prepare\_GWAS\_allele\_trait\_table\_height\_converted\_OSR.py
This script generates hypothetical offspring sex ratios of participants from their standing height, and simulate the sexes of the siblings of the participants according to their hypothetical offspring sex ratios. The number of participants and their siblings in the simulation are the same as the numbers in the UKB. The procedure is shown in **Figure 3a**

### 06.2\_GWAS\_height\_converted\_OSR\_logistic.py
The script performs GWAS on hypothetically generated siblings' sex to find potential SNPs underlying the hypothetical offspring sex ratio. The results are shown in **Figure 3b** of the paper.

### 07.1\_SLiM\_stabilizing.slim
The script performs SLiM simulations of sex ratio evolution in humans under stabilizing selection on population sex ratio of 0.5. The results are shown in **Figure 4** of the paper.

### 07.2\_SLiM\_direction.slim
The script conducts SLiM simulations of sex ratio evolution in humans under directional selection, starting with an initial population sex ratio of 0.5 and targeting an optimal of 0.524. The results are shown in **Figure 4** of the paper.

### 99\_GWAS\_false\_positive\_analysis.ipynb
The script verifies that treating each sibling birth event as an independent event in the GWAS is statistically unbiased.
