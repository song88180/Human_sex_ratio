import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import argparse
import numpy as np
import time
import pandas as pd
from pgenlib import PgenReader
import statsmodels.api as sm
import scipy.stats as stats

parser = argparse.ArgumentParser(description="Process some parameters.")

parser.add_argument("--chr", type=str, help="Chromosome number")
parser.add_argument("--start", type=int, help="start snp index (starting from 0)")
parser.add_argument("--end", type=int, help="end snp index (not included)")

args = parser.parse_args()
Chr_idx = args.chr
snp_start_idx = args.start
snp_end_idx = args.end

pfile_PATH = "/scratch/sigbio_project_root/sigbio_project19/siliangs/UKB_imp_filtered_SR_RAP"
GWAS_PATH = "/home/siliangs/GWAS_file"

phe_df_ori = pd.read_table(f'{GWAS_PATH}/GWAS_adj_extended.phe')
phe_df_ori = phe_df_ori[['IID','N_full_brothers','N_full_sisters']]
cov_df_ori = pd.read_table(f'{GWAS_PATH}/GWAS_adj_01.cov')
cov_df_ori = cov_df_ori[['IID','Genetic_sex','YoB']+[f'PC_{i}' for i in range(1,11)]]

phe_df_ori = pd.merge(left=phe_df_ori,right=cov_df_ori,on='IID')
phe_df_ori['N_brosis'] = phe_df_ori['N_full_brothers'] + phe_df_ori['N_full_sisters']
phe_df_ori['SR_brosis'] = phe_df_ori['N_full_brothers']/phe_df_ori['N_brosis']

df_dict = pd.read_table(f'{GWAS_PATH}/IID_new_old_dict.csv')
phe_df_ori = pd.merge(right=df_dict, left=phe_df_ori, right_on='OLD_IID', left_on='IID')
phe_df_ori['IID'] = phe_df_ori['RAP_IID'] 
phe_df_ori = phe_df_ori.drop(columns=['RAP_IID','OLD_IID'])


res_dict = {
    'P':[],
    't':[],
    'SE':[],
    'BETA':[]
}

print(f"chr{Chr_idx}:{snp_start_idx}-{snp_end_idx}")

var_df = pd.read_table(f'{pfile_PATH}/chr{Chr_idx}.pvar')
sam_df = pd.read_table(f'{pfile_PATH}/chr{Chr_idx}.psam')

print(f"varaince: {len(var_df)}, sample: {len(sam_df)}")

phe_df = phe_df_ori.copy()

target_IID = sam_df['IID']

phe_df = pd.merge(right=phe_df,left=target_IID,on='IID')

assert len(phe_df) == len(sam_df)

N_sample = len(sam_df)

# make sure IID in phe_df and .psam file have the same order.
assert (phe_df['IID'] != target_IID).sum() == 0

t1 = time.time()

N_sis_raw = phe_df['N_full_sisters'].values
N_bro_raw = phe_df['N_full_brothers'].values

Y_raw = phe_df['SR_brosis'].values

X_raw = phe_df[['Genetic_sex','YoB'] + [f'PC_{i}' for i in range(1,11)]].values

X_phe = []
Y_phe = []
filter_idx = []

for i in range(len(N_bro_raw)):
    if np.isnan(N_bro_raw[i]) or np.isnan(N_sis_raw[i]) or (N_bro_raw[i] + N_sis_raw[i] == 0):
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


for i in range(snp_start_idx, snp_end_idx):

    if i % 10 == 0:
        print(i)

    if (i - snp_start_idx) % 1000 == 0:
        out_buff = np.empty(shape=[min(1000, snp_end_idx - i), N_sample], dtype=np.int8)

        with PgenReader(
            f"{pfile_PATH}/chr{Chr_idx}.pgen".encode("utf-8"),
            raw_sample_ct = N_sample
        )  as reader:
            # variant_ct = reader.get_variant_ct()
            reader.read_range(
                variant_idx_start = i,
                variant_idx_end = min(1000+i, snp_end_idx),
                geno_int_out = out_buff
            )

    #reader.read(
    #    variant_idx = i,
    #    geno_int_out = out_buff
    #)

    X = np.concatenate((out_buff[(i - snp_start_idx) % 1000, filter_idx].reshape(-1,1), X_phe), axis=1)

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

out_df = pd.concat([var_df.iloc[snp_start_idx:snp_end_idx, :].reset_index(drop=True), res_df], axis=1)

out_df.to_csv(f'./WLS_brosis_cov_PC_for_all_SNPs_logistic_result/WLS_chr{Chr_idx}_{snp_start_idx}_{snp_end_idx}.tsv',sep='\t',index=False)
