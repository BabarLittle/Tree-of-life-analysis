# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:25:44 2023

@author: tangu
"""

#%%
import matplotlib as plt


#%%

#what I have to do before
mu=0.1
Q_matrix=np.array([[-3,1,1,1],[1,-3,1,1],[1,1,-3,1],[1,1,1,-3]])*mu
prob_dict=dict()
node_dict=dict()

#%%


probability=nucleotidic_score(msa_test1,test_1,Q_matrix,prob_dict,node_dict)

#to calculate the lnL sum per nucleotide
score=np.sum(probability["231"],axis=1)*0.25
np.log(sum(score))

#%%

#for the test_1

mu=[x/50.0 for x in range(1,21)]

scores=[]

for i in mu:
    Q_matrix=np.array([[-3,1,1,1],[1,-3,1,1],[1,1,-3,1],[1,1,1,-3]])*i
    prob_dict=dict()
    node_dict=dict()
    probability=nucleotidic_score(msa_test1,test_1,Q_matrix,prob_dict,node_dict)
    score=np.log(sum(np.sum(probability["231"],axis=1)*0.25))
    scores.append(score)
    
plt.plot(mu,scores)
plt.title("x/50")
plt.show()
print(f"The ideal mu is {mu[scores.index(max(scores))]}")
