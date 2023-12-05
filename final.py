# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:25:44 2023

@author: tangu
"""

#%%
import matplotlib.pyplot as plt
import numpy as np
import final_sequence_similarity as fss
import newick_form_main as nfm


#%%
branch_test1 = "data/ENSG00000013016_EHD3_NT.branchlength.dat"
msa_test1 = "data/ENSG00000013016_EHD3_NT.msa.dat"
table_test1 = "data/ENSG00000013016_EHD3_NT.table.dat"

branch_test2 = "data/ENSG00000112282_MED23_NT.branchlength.dat"
msa_test2 = "data/ENSG00000112282_MED23_NT.msa.dat"
table_test2 = "data/ENSG00000112282_MED23_NT.table.dat"

branch_test3 = "data/ENSG00000112984_KIF20A_NT.branchlength.dat"
msa_test3 = "data/ENSG00000112984_KIF20A_NT.msa.dat"
table_test3 = "data/ENSG00000112984_KIF20A_NT.table.dat"

test_1=nfm.newick_form(branch_test1,msa_test1,table_test1)
test_2=nfm.newick_form(branch_test2,msa_test2,table_test2)
test_3=nfm.newick_form(branch_test3,msa_test3,table_test3)
#%%

#what I have to do before
mu=0.1
Q_matrix=np.array([[-3,1,1,1],[1,-3,1,1],[1,1,-3,1],[1,1,1,-3]])*mu
prob_dict=dict()
node_dict=dict()

#%%


probability=fss.nucleotidic_score(msa_test1,test_1,Q_matrix,prob_dict,node_dict, mu)

#to calculate the lnL sum per nucleotide
score=np.sum(probability["231"],axis=1)*0.25
np.log(sum(score))

#%%

#for the test_1

mu=[x/50.0 for x in range(1,21)]



scores=[]
iteration = 0

for i in mu:
    print(f"iteration {iteration} / {len(mu)}")
    
    Q_matrix=np.array([[-3,1,1,1],[1,-3,1,1],[1,1,-3,1],[1,1,1,-3]])*i
    prob_dict=dict()
    node_dict=dict()
    probability=fss.nucleotidic_score(msa_test1,test_1,Q_matrix,prob_dict,node_dict, mu)
    score=np.log(sum(np.sum(probability["231"],axis=1)*0.25))
    scores.append(score)
    iteration += 1
print(f"iteration {iteration} / {len(mu)}")
    
plt.plot(mu,scores)
plt.title("x/50")
plt.show()
print(f"The ideal mu is {mu[scores.index(max(scores))]}")
