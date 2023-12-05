#Bonjour
import common as imp
import os
import numpy as np

branch_test = "data/ENSG00000013016_EHD3_NT.branchlength.dat"
msa_test = "data/ENSG00000013016_EHD3_NT.msa.dat"
table_test = "data/ENSG00000013016_EHD3_NT.table.dat"


tree = imp.construct_tree(branch_test, msa_test, table_test)
probability = imp.prob_calculation(tree)
likelihood = sum(sum(np.log(probability)))
print(likelihood)
