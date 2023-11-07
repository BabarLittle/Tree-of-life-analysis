# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 15:57:30 2023

@author: tangu
"""

def newick_form(branchlength, msa, table):
    bl=open(branchlength,"r")
    msa=open(msa,"r")
    table=open(table,"r")
    
    bl=bl.readlines()[0].split(",")
    msa=dict(line.replace("\n","").split(" ") for line in msa.readlines())
    table=[pos.replace("\n","").split(",") for pos in table.readlines()]
    
    newick={}
    
    #index for branch length and sequence
    index=0
    
    for name in table:
        
        #name [0] is the parent, name[1] is the child
        
        #for the parent name[0]
        if name[0] not in newick.keys(): 
            #create a dictionary as values of name (name is keys of newick dict)
            newick[name[0]]={}
            #create an empty list for children
            newick[name[0]]["children"]=[]
                
        
        #add the child name[1] of name[0] in the list
        newick[name[0]]["children"].append(name[1])
          
        #only conditions needed for name[1] is to know if it is already in newick.keys()
        #to create an empty dict and to add sequence (because only last individual have sequence), because it would be only once children
        #so we have to add the only parent, branch and length every time
        if name[1] not in newick.keys():
            #create a dictionary as values of name[1]
            newick[name[1]]={}
            #add the sequence of i[1]
            newick[name[1]]["sequence"]=msa[name[1]]
        
        #add the parent of i[1] (not in a list because only 1 parent)
        newick[name[1]]["parent"]=name[0]
        #add the value of the branch length
        newick[name[1]]["branch_length"]=bl[index]

        
        
        index+=1
        
    
    return(newick)
    


