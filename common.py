import os

def construct_tree(branchlength, msa, table ):
    #initialize the dictionary
    tree = {}
    list_potential_key = [] 
    
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
                    
                    index += 1          
                 
    except IOError:
        print("Error file not found")
    return tree



