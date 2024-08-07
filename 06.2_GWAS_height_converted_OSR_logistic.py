\import numpy as np
import pandas as pd
from pgenlib import PgenReader
import statsmodels.api as sm
import time
import sys

# Get the replication number from command-line arguments
rep = int(sys.argv[1])

# Load the SNP table for the given replication
df_all = pd.read_table(f'./SNP_table_OSR_height_resid/SNP_table_OSR_height_resid_rep{rep}.tsv')

# Initialize dictionaries to store data
N_list = {}
y_list = {}
N_bro = {}
N_sis = {}

# Process data for each allele count (0, 1, 2)
for i in range(3):

    N_list[i] = list(df_all[f'{i}_bro']+df_all[f'{i}_sis'])
    N_bro[i] = df_all[f'{i}_bro']
    N_sis[i] = df_all[f'{i}_sis']
    y_list[i] = (N_bro[i]/np.array(N_list[i])).fillna(0) # Handle cases where N == 0

# Initialize lists to store regression results
beta_list = []
se_list = []
P_list = []
t_list = []

# Create design matrix for regression
X = np.array([0,1,2]).reshape(-1,1)
X = sm.add_constant(X)

# Perform logistic regression for each SNP
for i in range(len(df_all)):
    
    if i % 10000 == 0:
        print(i,end='\r',flush=True)
    
    # Logistic regression using proportion and weighting
    glm = sm.GLM(
        [y_list[0][i],y_list[1][i],y_list[2][i]],
        X,
        family=sm.families.Binomial(),
        freq_weights=[N_list[0][i], N_list[1][i], N_list[2][i]]
    )

    res = glm.fit()
    
    beta_list.append(res.params[1])
    se_list.append(res.bse[1])
    P_list.append(res.pvalues[1])
    t_list.append(res.tvalues[1])

df_all['P'] = P_list
df_all['t'] = t_list
df_all['SE'] = se_list
df_all['BETA'] = beta_list

# Save the results
df_all.to_csv(f'./GWAS_OSR_from_height_resid/GWAS_OSR_from_height_resid_rep{rep}.tsv', index=False, sep='\t')
