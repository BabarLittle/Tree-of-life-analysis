import os
import math
import numpy as np


branch_test = "data\ENSG00000013016_EHD3_NT.branchlength.dat"
msa_test = "data/ENSG00000013016_EHD3_NT.msa.dat"
table_test = "data/ENSG00000013016_EHD3_NT.table.dat"
mu = -0.1
Q_matrix=np.array([[-3,1,1,1],[1,-3,1,1],[1,1,-3,1],[1,1,1,-3]])*mu

def construct_tree(branchlength, msa, table ):
    #initialize the dictionary
    tree = {}
    
    #open file
    try:
        with open(branchlength, 'r') as lb, \
             open(msa, 'r') as msa, \
             open(table, 'r') as table:
                 
                 #store data in list format
                 len_branch = lb.readlines()[0].split(",")
                 msa_data = dict(line.split(" ") for line in msa)
                 
                 index = 0 #index to iterate over length branch list
                 for line in table:
                    
                    #Children part since the first part of the line is the parent
                    if line.split(",")[0] not in tree.keys(): #key check and creation 
                        tree[line.split(",")[0]] = {}
                        #adding children
                        tree[line.split(",")[0]]["children"] = [line.split(",")[1].strip("\n")]
                    else:
                        #adding children with existing parent key
                        tree[line.split(",")[0]]["children"].append(line.split(",")[1].strip("\n"))
                    
                    # Parent related part since the last part of the line is the children
                    if line.split(",")[1].strip("\n") not in tree.keys(): #terminated node
                
                        tree[line.split(",")[1].strip("\n")] = {}
                        tree[line.split(",")[1].strip("\n")]["Parents"] = line.split(",")[0] #set parents
                        tree[line.split(",")[1].strip("\n")]["Length_to_P"] = len_branch[index] # set length to parents
                        tree[line.split(",")[1].strip("\n")]["Sequence"] = msa_data[line.split(",")[1].strip("\n")] # set sequence
                        
                    else: #is related to a parent
                        tree[line.split(",")[1].strip("\n")]["Parents"] = line.split(",")[0] #set parents
                        tree[line.split(",")[1].strip("\n")]["Length_to_P"] = len_branch[index] #set length
                    
                    index += 1 # increase index      
                 
    except IOError:
        print("Error file not found")
    return tree



def alignment_score(sequence, branch_length):
    result_vector = []
    for i in range(len(sequence)):
                    if sequence[i]=="A": vector = [1,0,0,0]
                    elif sequence[i]=="C": vector = [0,1,0,0]
                    elif sequence[i]=="G": vector = [0,0,1,0]
                    elif sequence[i]=="T": vector = [0,0,0,1]
                    
                    prob = np.matmul(np.exp(Q_matrix*branch_length), vector)
                    result_vector.append(prob)
    return math.prod(result_vector) 

def alignment_similarity(species1, species2):
    return math.prod([species1, species2])


def prob_calculation(dict_tree, dynamic_dict = None, keys = None):
    #Initialize the dictionary and list if not given
    if dynamic_dict is None:
        dynamic_dict = {}
    if keys is None:
        keys = list(dict_tree.keys())
    
     # Get the first element containing "_" in the list of keys
    try:
        index = next(i for i, e in enumerate(keys) if "_" in e)
        sequence = dict_tree[keys[index]]["Sequence"] #Get sequence
        if "Length_to_P" in dict_tree[keys[index]]: #check if length to parent is in the nested dictionary
            branch_length = float(dict_tree[keys[index]]["Length_to_P"]) #transform length to float
            
            dynamic_dict[keys[index]] = alignment_score(sequence, branch_length) #calculate alignment score
            keys.remove(keys[index]) #remove the key from the list of keys

            
        return prob_calculation(dict_tree, dynamic_dict=dynamic_dict, keys = keys) #recursive call
    
        
    except StopIteration:
        pass
    
    #calculate alignment similarity of parent based on children
    for key in keys: #for each key in the list of keys
        children = tree[key]["children"] #get children
        if all(child in dynamic_dict for child in children): #check if children are in the dynamic dictionary
            dynamic_dict[key] = alignment_similarity(dynamic_dict[children[0]], dynamic_dict[children[1]]) #calculate alignment similarity
            keys.remove(key) #remove the key from the list of keys
     
    #if there are still keys in the list of keys, recursive call   
    if len(keys) > 0:
        return prob_calculation(dict_tree, dynamic_dict=dynamic_dict, keys = keys) #recursive call
    else: #return result
        return dynamic_dict
        

        
        
            
    
        
        
        
        
tree = construct_tree(branch_test, msa_test, table_test)
likelyhood = sum(prob_calculation(tree)['231'])
print(likelyhood)




    



