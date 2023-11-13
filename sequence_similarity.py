# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:22:48 2023

@author: tangu
"""

#%%
import math
import numpy as np
import random



#%%

#create alignment probability for each nucleotide (only for species with fasta sequence)

def alignment_probability(species,mu,test):
    
    Q_matrix=np.array([[-3,1,1,1],[1,-3,1,1],[1,1,-3,1],[1,1,1,-3]])*mu
    seq=test[species]["sequence"]
    bl=test[species]["branch_length"]
    prob_alignment=[]
    
    for i in range(len(seq)):
        if seq[i]=="A": vec=[1,0,0,0]
        elif seq[i]=="C": vec=[0,1,0,0]
        elif seq[i]=="G": vec=[0,0,1,0]
        elif seq[i]=="T": vec=[0,0,0,1]
        
        prob=np.matmul(np.exp(Q_matrix*bl),vec)
        prob_alignment.append(list(prob))
    
    prob_dict[species]=prob_alignment
    

#%%

#create alignment probability of the parent of 2 individuals
#based on their alignment probability


def alignment_similarity(spe_1,spe_2,prob_dict,test):
    if test[spe_1]["parent"]==test[spe_2]["parent"]:
        prob_ancester=[]
        prob_1=prob_dict[spe_1]
        prob_2=prob_dict[spe_2]
        for j in range(len(prob_1)):
            prob_ancester.append([prob_1[j][i]*prob_2[j][i] for i in range(4)])
        prob_dict[test[spe_1]["parent"]]=prob_ancester




#%%  
#open msa
 
def msa_opening(msa_test):
    msa=open(msa_test,"r")
    msa=dict(line.replace("\n","").split(" ") for line in msa.readlines())
    
    return(msa)

#%%

#create dictionnary with all alignment probability from sequence we have

def probability_dict(msa,test_1,mu,prob_dict):
    names=list(msa.keys())
    for i in names:
        alignment_probability(i, mu, test_1)
    
    return(prob_dict) 
    
    
#%%

#choose randomly the first species to start with
  
def start_node(msa):
    ch=random.choice(list(msa.keys()))
    return(ch)


#%%

#find the single parent of an individual and then change parent variable

def find_parent(test,ch):
    parent=test[ch]["parent"]
    node_dict[parent]=test[parent]["children"]    
    
    return(parent)
    



#%%

#find both children of an individual

def find_children(test,parent):
    node_dict[parent]=test[parent]["children"]    

#%%

#check if teh alignment probability of a parent can be found
#can be found if alignment probability of both children are already found

def probability_vec(node_dict,prob_dict,parent):
    num=[0,0]
    for i in node_dict[parent]:
        if i in list(prob_dict.keys()):
            num[node_dict[parent].index(i)]=1
    
    return(num)    

