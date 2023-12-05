# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 01:36:33 2023

@author: tangu
"""

from collections import defaultdict
from treelib import Node, Tree

#without branch length


roots="226"
#231 or 229 or 226

# Create tree and add root node
tree = Tree()
tree.create_node(roots, roots)


# Recursive function to add children to tree
def add_children(parent):
    for child in test_3[parent]["children"]:
        try:
            tree.create_node(child, child, parent=parent)
            add_children(child)
        except KeyError:
            continue

# Add children to tree
add_children(roots)

tree.save2file("tree_test_3.txt")

