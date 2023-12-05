# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:22:48 2023

@author: tangu
@corrected by: jikael
"""

#%%
import math
import numpy as np
import random
from scipy.linalg import expm





#%%

#create alignment probability for each nucleotide (only for species with fasta sequence)

def alignment_probability(species,mu,test,prob_dict,Q_matrix):
    
    seq=test[species]["sequence"]
    bl=test[species]["branch_length"]
    prob_alignment=[]
    
    for i in range(len(seq)):
        if seq[i]=="A": vec=[1,0,0,0]
        elif seq[i]=="C": vec=[0,1,0,0]
        elif seq[i]=="G": vec=[0,0,1,0]
        elif seq[i]=="T": vec=[0,0,0,1]
        
        prob=list(vec)
        prob_alignment.append(prob)
    
    prob_dict[species]=prob_alignment
    

#%%

#create alignment probability of the parent of 2 individuals
#based on their alignment probability


def alignment_similarity(spe_1,spe_2,prob_dict,test,Q_matrix):
    if test[spe_1]["parent"]==test[spe_2]["parent"]:
        prob_ancester=[]
        prob_1=prob_dict[spe_1][:]
        prob_2=prob_dict[spe_2][:]
        bl_1=float(test[spe_1]["branch_length"])
        bl_2=float(test[spe_2]["branch_length"])

        
        for j in range(len(prob_1)):
            
            vec1=list(np.matmul(expm(Q_matrix*bl_1),prob_1[j]))
            vec2=list(np.matmul(expm(Q_matrix*bl_2),prob_2[j]))
            prob_ancester.append([vec1[i]*vec2[i] for i in range(4)])
            #is it really a multiplication (then big value for each ancesters)
            
            #this line is to scale the score in order to compare doesn't matter which generation it is
            #prob_ancester[j]=[prob_ancester[j][i]/sum(prob_ancester[j]) for i in range(4)]
        
        prob_dict[test[spe_1]["parent"]]=prob_ancester




#%%  
#open msa
 
def msa_opening(msa_test):
    msa=open(msa_test,"r")
    msa=dict(line.replace("\n","").split(" ") for line in msa.readlines())
    
    return(msa)

#%%

#create dictionnary with all alignment probability from sequence we have

def probability_dict(msa,test_1,mu,prob_dict,Q_matrix):
    names=list(msa.keys())
    for i in names:
        alignment_probability(i, mu, test_1,prob_dict,Q_matrix)
    
    return(prob_dict) 
    
    
#%%

#choose randomly the first species to start with
  
def start_node(msa):
    ch=random.choice(list(msa.keys()))
    return(ch)


#%%

#find the single parent of an individual and then change parent variable

def find_parent(test,ch, node_dict):
    parent=test[ch]["parent"]
    node_dict[parent]=test[parent]["children"]    
    
    return(parent)
    



#%%

#find both children of an individual

def find_children(test,parent, node_dict):
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

#%%

def nucleotidic_score(msa_test,test,Q_matrix,prob_dict,node_dict, mu):

    #open msa
    msa=msa_opening(msa_test) 
        
    #create dictionnary with all alignment probability from sequence we have
    prob_dict=probability_dict(msa,test,mu,prob_dict,Q_matrix)
    #choose randomly the first species to start with
    start=start_node(msa)
        
    #get the parent and family from the starter
    parent=find_parent(test,start, node_dict)
    #loop that will start here
        
    for i in range(1000):
        #vector that check if the children's probability alignment are in the prob_dict
        num=probability_vec(node_dict,prob_dict,parent)
            
        #create alignment probability if possible or get down of one branch
        species=node_dict[parent]
        #if num=[1,1] then both alignments are in prob dict, it is possible to get
        #the alignment of the ancester
        if num==[1,1]:
            alignment_similarity(species[0],species[1],prob_dict,test,Q_matrix)
            del node_dict[parent]
            try:
                parent=find_parent(test,parent, node_dict)
            except KeyError:
                break
        #if num!=[1,1], it means that at least one of the alignment is not in prob dict
        #so we have to go further in the tree to find it.
        else:
            parent=species[num.index(0)] #parent become the 0 in probability vector
            find_children(test,parent, node_dict) #add the family of the new parent
    return(prob_dict)