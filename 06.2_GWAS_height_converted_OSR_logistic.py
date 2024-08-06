import numpy as np
import time
import pandas as pd
from pgenlib import PgenReader
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
import statsmodels.formula.api as smf
import scipy.stats as stats
import sys

rep = int(sys.argv[1])

df_all = pd.read_table(f'../aggregation_method/Allele_table_from_height_resid/Allele_table_from_height_resid_rep{rep}.tsv')

N_list = {}
y_list = {}
N_bro = {}
N_sis = {}
for i in range(3):

    N_list[i] = list(df_all[f'{i}_bro']+df_all[f'{i}_sis'])
    N_bro[i] = df_all[f'{i}_bro']
    N_sis[i] = df_all[f'{i}_sis']
    y_list[i] = (N_bro[i]/np.array(N_list[i])).fillna(0) # N == 0

reg = LinearRegression()

beta_list = []
se_list = []
P_list = []
t_list = []


X = np.array([0,1,2]).reshape(-1,1)
X = sm.add_constant(X)


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

df_all.to_csv(f'./Sex_ratio_from_height_resid/WLS_brosis_from_height_resid_logistic_rep{rep}.tsv', index=False, sep='\t')
