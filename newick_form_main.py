# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 18:21:31 2023

@author: tangu
"""

branch_test1 = "ENSG00000013016_EHD3_NT.branchlength.dat"
msa_test1 = "ENSG00000013016_EHD3_NT.msa.dat"
table_test1 = "ENSG00000013016_EHD3_NT.table.dat"

branch_test2 = "ENSG00000112282_MED23_NT.branchlength.dat"
msa_test2 = "ENSG00000112282_MED23_NT.msa.dat"
table_test2 = "ENSG00000112282_MED23_NT.table.dat"

branch_test3 = "ENSG00000112984_KIF20A_NT.branchlength.dat"
msa_test3 = "ENSG00000112984_KIF20A_NT.msa.dat"
table_test3 = "ENSG00000112984_KIF20A_NT.table.dat"


test_1=newick_form(branch_test1,msa_test1,table_test1)
test_2=newick_form(branch_test2,msa_test2,table_test2)
test_3=newick_form(branch_test3,msa_test3,table_test3)