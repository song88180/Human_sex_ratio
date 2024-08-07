import numpy as np
import numpy.random as nrand
import time
import pandas as pd
from pgenlib import PgenReader
import statsmodels.api as sm 
import sys

# Get the replication number from command-line arguments
rep = int(sys.argv[1])

# Load and preprocess phenotype and covariate data
phe_df_raw = pd.read_table('/home/siliang/Public/workspace/4_UKB_merged_files/GWAS/GWAS_adj_extended.phe')
phe_df_raw = phe_df_raw[['FID','IID','Sex_merged','N_full_brothers','N_full_sisters','Standing_height']]
phe_df_raw['SR_from_height'] = phe_df_raw['Standing_height']/(2*phe_df_raw['Standing_height'].mean())
cov_df = pd.read_table('/home/siliang/Public/workspace/4_UKB_merged_files/GWAS/GWAS_adj_01.cov')

# Merge phenotype and covariate data
phe_df_raw = pd.merge(left = phe_df_raw, right = cov_df[['IID','Genetic_sex','Age','Age_squared']+[f'PC_{i}' for i in range(1,11)]], on='IID', copy=False)

# Filter out rows with NA values
idx_notna = phe_df_raw[['Genetic_sex', 'Age', 'Age_squared', 'SR_from_height']+[f'PC_{i}' for i in range(1,11)]].isna().sum(axis=1) == 0

# Regress out covariates from SR_from_height
X = phe_df_raw.loc[idx_notna, ['Genetic_sex', 'Age', 'Age_squared']+[f'PC_{i}' for i in range(1,11)]]
y = phe_df_raw.loc[idx_notna, 'SR_from_height'] 
X = sm.add_constant(X) 
est = sm.OLS(y, X).fit() 
phe_df_raw.loc[idx_notna, 'SR_from_height_resid']  = est.resid + 0.5

# Set NA values to 0 for specified columns
_idx = phe_df_raw[['N_full_brothers','N_full_sisters','SR_from_height_resid']].isna().sum(axis=1) > 0
phe_df_raw.loc[_idx,['N_full_brothers','N_full_sisters','SR_from_height_resid']] = 0.0

# Simulate number of brothers and sisters based on height residuals
N_brothers_from_height = nrand.binomial(
    n=(phe_df_raw['N_full_brothers']+phe_df_raw['N_full_sisters']),
    p=phe_df_raw['SR_from_height_resid']
)
N_sisters_from_height = phe_df_raw['N_full_brothers']+phe_df_raw['N_full_sisters'] - N_brothers_from_height

phe_df_raw['N_brothers_from_height'] = N_brothers_from_height
phe_df_raw['N_sisters_from_height'] = N_sisters_from_height

# Process each chromosome
for Chr_idx in list(range(1,23))+['X']:

    print(f"chr{Chr_idx}")
    
    # Load variant and sample data
    var_df = pd.read_table(f'/home/siliang/Public/workspace/UKB_genome_tmp/chr{Chr_idx}/chr{Chr_idx}.pvar')
    sam_df = pd.read_table(f'/home/siliang/Public/workspace/UKB_genome_tmp/chr{Chr_idx}/chr{Chr_idx}.psam')
    
    print(f"varaince: {len(var_df)}, sample: {len(sam_df)}")
    
    target_IID = sam_df['IID']
    
    # Merge sample data with phenotype data
    phe_df = pd.merge(right=phe_df_raw,left=target_IID,on='IID')
    
    # Ensure IID in phe_df and sam_df have the same order
    assert len(phe_df) == len(sam_df)

    # Ensure IID in phe_df and .psam file have the same order.
    assert (phe_df['IID'] != target_IID).sum() == 0

    N_sample = len(sam_df)
    
    t1 = time.time()

    # Prepare matrix B for faster computation
    B = phe_df.loc[:,['N_brothers_from_height','N_sisters_from_height']].to_numpy()
    
    bro_dict = {0: [], 1: [], 2: []}
    sis_dict = {0: [], 1: [], 2: []}

    with PgenReader(
        f"/home/siliang/Public/workspace/UKB_genome_tmp/chr{Chr_idx}/chr{Chr_idx}.pgen".encode("utf-8"),
        raw_sample_ct=N_sample
    )  as reader:
        variant_ct = reader.get_variant_ct()

        # initialize the output buffer
        out_buff = np.empty(N_sample, dtype=np.int8)

        for i in range(variant_ct):

            if i % 10000 == 0:
                print(i, end='\r', flush=True)

            reader.read(
                variant_idx=i,
                geno_int_out=out_buff
            )

            # By default, 0/1/2 represents the alternate allele count.
            for N_allele in [0,1,2]:
                A = (out_buff == N_allele)
                C = A.dot(B)   # This step calculates total number of brothers/sisters of participants with genotype 0/1/2.
                bro_dict[N_allele].append(C[0])
                sis_dict[N_allele].append(C[1])
                if sum(C) == 0:
                    continue

    t2 = time.time()
    
    print(f"total time: {time.strftime('%H:%M:%S', time.gmtime(t2-t1))}")

    var_df['0_bro'] = bro_dict[0]
    var_df['0_sis'] = sis_dict[0]
    var_df['1_bro'] = bro_dict[1]
    var_df['1_sis'] = sis_dict[1]
    var_df['2_bro'] = bro_dict[2]
    var_df['2_sis'] = sis_dict[2]

    var_df = var_df.astype(
        {
            '0_bro': int,
            '0_sis': int,
            '1_bro': int,
            '1_sis': int,
            '2_bro': int,
            '2_sis': int
        }
    )
    
    # Save results
    var_df.to_csv(f'SNP_table_OSR_height_resid/SNP_table_OSR_height_resid_chr{Chr_idx}_rep{rep}.tsv',sep='\t',index=False)
