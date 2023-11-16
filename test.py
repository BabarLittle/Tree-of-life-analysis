# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 13:25:44 2023

@author: tangu
"""


#%%

#what I have to do before
mu=-0
Q_matrix=np.array([[-3,1,1,1],[1,-3,1,1],[1,1,-3,1],[1,1,1,-3]])*mu
prob_dict=dict()
node_dict=dict()
#open msa
msa=msa_opening(msa_test2) 
    
#create dictionnary with all alignment probability from sequence we have
prob_dict=probability_dict(msa,test_2,mu,prob_dict,Q_matrix)
print(prob_dict.keys())
print(len(prob_dict.keys()))
#choose randomly the first species to start with
start=start_node(msa)
print(start)
    
#get the parent and family from the starter
parent=find_parent(test_2,start)
print(parent)
print(node_dict)
#loop that will start here
    
for i in range(1000):
    #vector that check if the children's probability alignment are in the prob_dict
    num=probability_vec(node_dict,prob_dict,parent)
    print(num)
        
    #create alignment probability if possible or get down of one branch
    species=node_dict[parent]
    print(species)
    #if num=[1,1] then both alignments are in prob dict, it is possible to get
    #the alignment of the ancester
    if num==[1,1]:
        alignment_similarity(species[0],species[1],prob_dict,test_2,Q_matrix)
        del node_dict[parent]
        try:
            parent=find_parent(test_2,parent)
        except KeyError:
            break
        print(node_dict)
    #if num!=[1,1], it means that at least one of the alignment is not in prob dict
    #so we have to go further in the tree to find it.
    else:
        parent=species[num.index(0)] #parent become the 0 in probability vector
        print(parent)
        find_children(test_2,parent) #add the family of the new parent
        print(node_dict)

print(prob_dict["Erinaceus_europaeus"])
print(test_2["Erinaceus_europaeus"]["parent"])
print(prob_dict["181"])
print(prob_dict["229"])
#print(prob_dict.keys())
