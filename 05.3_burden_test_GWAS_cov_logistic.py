import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import argparse
import time
import numpy as np
import pandas as pd
from pgenlib import PgenReader
import statsmodels.api as sm

parser = argparse.ArgumentParser(description="Process some parameters.")

parser.add_argument("--chr", type=str, help="Chromosome number")
args = parser.parse_args()
Chr_idx = args.chr

res_dict = {
    'P':[],
    't':[],
    'SE':[],
    'BETA':[]
}

phe_df_raw = pd.read_table('/home/siliang/Public/workspace/4_UKB_merged_files/GWAS/GWAS_adj_extended_new.phe')
phe_df_raw = phe_df_raw[['IID','N_full_brothers','N_full_sisters','Sex_merged']]
_idx = phe_df_raw[['N_full_brothers','N_full_sisters']].isna().sum(axis=1) > 0
phe_df_raw.loc[_idx,['N_full_brothers','N_full_sisters']] = 0.0

cov_df_raw = pd.read_table('/home/siliang/Public/workspace/4_UKB_merged_files/GWAS/GWAS_adj_01_new.cov')
cov_df_raw = cov_df_raw[['IID','Genetic_sex','YoB']+[f'PC_{i}' for i in range(1,11)]]

phe_df_raw = pd.merge(left=phe_df_raw,right=cov_df_raw,on='IID')

sam_KM_df = pd.read_table('/home/siliang/Public/workspace/UKB_imp_filtered_SR_RAP/chr1.psam')

print(f"chr{Chr_idx}")

var_WES_df = pd.read_table(f'./mask_bfiles_0.05_0.01_0.001/step2_chr{Chr_idx}_masks.bim',names=['chr','name','0','pos','class','ref'])
sam_WES_df = pd.read_table(f'./mask_bfiles_0.05_0.01_0.001/step2_chr{Chr_idx}_masks.fam',names=['FID','IID','0','1','sex','-9'])

print(f"variant: {len(var_WES_df)}, sample: {len(sam_WES_df)}")

sam_KM_WES_df = pd.merge(left=sam_WES_df['IID'], right=sam_KM_df, on='IID', how='left')
# make sure IID in sam_KM_WES_df and in sam_WES_df have the same order.
assert (sam_KM_WES_df['IID'] != sam_WES_df['IID']).sum() == 0

# creat mask to only include European ancestry individuals, and exclude individuals whose siblings have appeared in the analysis.
KM_WES_mask = np.array(sam_KM_WES_df['#FID'].notna())

KM_WES_phe_df = pd.merge(left=sam_WES_df['IID'], right=phe_df_raw, on='IID', how='left')

# make sure IID in KM_WES_phe_df and in sam_WES_df have the same order.
assert (KM_WES_phe_df['IID'] != sam_WES_df['IID']).sum() == 0

N_sample = len(sam_WES_df)

t1 = time.time()

X_raw = KM_WES_phe_df[['Genetic_sex','YoB'] + [f'PC_{i}' for i in range(1,11)]].values
X_phe = []
Y_phe = []
filter_idx = []

N_sis_raw = KM_WES_phe_df['N_full_sisters'].values
N_bro_raw = KM_WES_phe_df['N_full_brothers'].values

for i in range(len(N_bro_raw)):
    if np.isnan(N_bro_raw[i]) or np.isnan(N_sis_raw[i]) or (N_bro_raw[i] + N_sis_raw[i] == 0) or (not KM_WES_mask[i]):
        continue
    for j in range(int(N_bro_raw[i])):
        X_phe.append(X_raw[i,:].tolist())
        Y_phe.append(1)
        filter_idx.append(i)
    for j in range(int(N_sis_raw[i])):
        X_phe.append(X_raw[i,:].tolist())
        Y_phe.append(0)
        filter_idx.append(i)
Y_phe = np.array(Y_phe)
X_phe = np.array(X_phe)

# initialize the output buffer
out_buff = np.empty(N_sample, dtype=np.int8)

variant_ct = len(var_WES_df)

for i in range(variant_ct):

    if i % 10 == 0:
        print(i,end='\r',flush=True)

    if var_WES_df['class'].iloc[i] != 'M3.0.01':
        continue

    with PgenReader(
        f"./mask_bfiles_0.05_0.01_0.001/step2_chr{Chr_idx}_masks.bed".encode("utf-8"),
        raw_sample_ct = N_sample
    )  as reader:

        reader.read(
            variant_idx = i,
            geno_int_out = out_buff
        )

    X = np.concatenate((out_buff[filter_idx].reshape(-1,1), X_phe), axis=1) 

    idx_ = (X[:,0] != -9)

    Y = Y_phe[idx_]

    X = X[idx_]

    X = sm.add_constant(X)

    try:
        mod = sm.Logit(Y, X)
        res = mod.fit(disp=False)

        res_dict['BETA'].append(res.params[1])
        res_dict['SE'].append(res.bse[1])
        res_dict['P'].append(res.pvalues[1])
        res_dict['t'].append(res.tvalues[1])

    except Exception as e:
        res_dict['BETA'].append(np.nan)
        res_dict['SE'].append(np.nan)
        res_dict['P'].append(np.nan)
        res_dict['t'].append(np.nan)
        continue

t2 = time.time()
print(f"total time: {time.strftime('%H:%M:%S', time.gmtime(t2-t1))}")

res_df = pd.DataFrame(res_dict)

out_df = pd.concat([var_WES_df.loc[var_WES_df['class'] == 'M3.0.01'].reset_index(drop=True), res_df], axis=1)

out_df.to_csv(f'./Gene_brosis_cov_logistic/Gene_brosis_cov_logistic_chr{Chr_idx}.tsv',sep='\t',index=False)
