import numpy as np
import numpy.random as nrand
import time
import pandas as pd
from pgenlib import PgenReader
import statsmodels.api as sm 
import sys
import os

rep = int(sys.argv[1])

phe_df_raw = pd.read_table('/home/siliang/Public/workspace/4_UKB_merged_files/GWAS/GWAS_adj_extended.phe')
phe_df_raw = phe_df_raw[['FID','IID','Sex_merged','N_full_brothers','N_full_sisters','Standing_height']]
phe_df_raw['SR_from_height'] = phe_df_raw['Standing_height']/(2*phe_df_raw['Standing_height'].mean())

cov_df = pd.read_table('/home/siliang/Public/workspace/4_UKB_merged_files/GWAS/GWAS_adj_01.cov')

phe_df_raw = pd.merge(left = phe_df_raw, right = cov_df[['IID','Genetic_sex','Age','Age_squared']+[f'PC_{i}' for i in range(1,11)]], on='IID', copy=False)

idx_notna = phe_df_raw[['Genetic_sex', 'Age', 'Age_squared', 'SR_from_height']+[f'PC_{i}' for i in range(1,11)]].isna().sum(axis=1) == 0

X = phe_df_raw.loc[idx_notna, ['Genetic_sex', 'Age', 'Age_squared']+[f'PC_{i}' for i in range(1,11)]]
y = phe_df_raw.loc[idx_notna, 'SR_from_height'] 
X = sm.add_constant(X) 
est = sm.OLS(y, X).fit() 
phe_df_raw.loc[idx_notna, 'SR_from_height_resid']  = est.resid + 0.5


_idx = phe_df_raw[['N_full_brothers','N_full_sisters','SR_from_height_resid']].isna().sum(axis=1) > 0
phe_df_raw.loc[_idx,['N_full_brothers','N_full_sisters','SR_from_height_resid']] = 0.0

N_brothers_from_height = nrand.binomial(
    n=(phe_df_raw['N_full_brothers']+phe_df_raw['N_full_sisters']),
    p=phe_df_raw['SR_from_height_resid'])

N_sisters_from_height = phe_df_raw['N_full_brothers']+phe_df_raw['N_full_sisters'] - N_brothers_from_height

phe_df_raw['N_brothers_from_height'] = N_brothers_from_height
phe_df_raw['N_sisters_from_height'] = N_sisters_from_height


for Chr_idx in list(range(1,23))+['X']:

    print(f"chr{Chr_idx}")

    var_df = pd.read_table(f'/home/siliang/Public/workspace/UKB_genome_tmp/chr{Chr_idx}/chr{Chr_idx}.pvar')
    sam_df = pd.read_table(f'/home/siliang/Public/workspace/UKB_genome_tmp/chr{Chr_idx}/chr{Chr_idx}.psam')

    print(f"varaince: {len(var_df)}, sample: {len(sam_df)}")

    target_IID = sam_df['IID']

    phe_df = pd.merge(right=phe_df_raw,left=target_IID,on='IID')

    # if either "N_full_brothers" or "N_full_sisters" is NaN, then assign 0 for both "N_full_brothers" and "N_full_sisters"
    #_idx = phe_df[['N_males_in_family','N_females_in_family']].isna().sum(axis=1) > 0
    #phe_df.loc[_idx,['N_males_in_family','N_females_in_family']] = 0.0

    assert len(phe_df) == len(sam_df)

    N_sample = len(sam_df)

    # make sure IID in phe_df and .psam file have the same order.
    assert (phe_df['IID'] != target_IID).sum() == 0

    t1 = time.time()

    # For faster speed
    B = phe_df.loc[:,['N_brothers_from_height','N_sisters_from_height']].to_numpy()

    P_list = []
    Chi2_list = []
    bro_dict = {0:[],1:[],2:[]}
    sis_dict = {0:[],1:[],2:[]}

    with PgenReader(
        f"/home/siliang/Public/workspace/UKB_genome_tmp/chr{Chr_idx}/chr{Chr_idx}.pgen".encode("utf-8"),
        raw_sample_ct = N_sample
    )  as reader:
        #a = reader.get_raw_sample_ct()
        variant_ct = reader.get_variant_ct()

        # initialize the output buffer
        out_buff = np.empty(N_sample, dtype=np.int8)

        for i in range(variant_ct):

            if i % 10000 == 0:
                print(i,end='\r',flush=True)

            reader.read(
                variant_idx = i,
                geno_int_out = out_buff
            )

            # by default, 0/1/2 represents the alternate allele count.
            for N_allele in [0,1,2]:
                A = (out_buff == N_allele)
                C = A.dot(B)
                bro_dict[N_allele].append(C[0])
                sis_dict[N_allele].append(C[1])
                if sum(C) == 0:
                    continue

            #res = chi2_contingency(N_bro_sis)
            #P_list.append(res.pvalue)
            #Chi2_list.append(res.statistic)

    t2 = time.time()

    #print(f"total time: {t2-t1}")
    print(f"total time: {time.strftime('%H:%M:%S', time.gmtime(t2-t1))}")

    var_df['0_bro'] = bro_dict[0]
    var_df['0_sis'] = sis_dict[0]
    var_df['1_bro'] = bro_dict[1]
    var_df['1_sis'] = sis_dict[1]
    var_df['2_bro'] = bro_dict[2]
    var_df['2_sis'] = sis_dict[2]

    var_df = var_df.astype(
        {
            '0_bro':int,
            '0_sis':int,
            '1_bro':int,
            '1_sis':int,
            '2_bro':int,
            '2_sis':int
        }
    )

    var_df.to_csv(f'Allele_table_from_height_resid/Allele_table_from_height_resid_chr{Chr_idx}_rep{rep}.tsv',sep='\t',index=False)


os.system(f'/bin/bash ./merge_Allele_table_from_height_resid.sh {rep}')
