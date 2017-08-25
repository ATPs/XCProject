# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 11:34:17 2015

@author: k
"""

class NcbiTaxonTree(object):
    """
    NCBI taxonomy tree class. 
    """
    def __init__(self, taxonid, _parent = None):
        self.taxonid = taxonid
        self._children = []
        self._parent = _parent
    
    def add_children(self, obj):
        self._children.append(obj)
    
#    def add_chldrens(self, obj):
#        
    
    def children(self):
        """
        return list of children
        """
        return self._children
    
    def childrenids(self):
        """
        return a list of children taxon ids
        """
        ids = []
        for child in self._children:
            ids.append(child.taxonid)
        return ids
    
    def parentid(self):
        """
        return parent id
        """
        if self._parent == None:
            return None
        return self._parent.taxonid
        
    def parent(self):
        """
        return parent
        """
        return self._parent
    
    def __str__(self):
        """
        return a string with the description of NCBI
        """
        return str(self.taxonid)


mylist = open("list2.txt","r").readlines()
for num in range(len(mylist)):
    mylist[num] = mylist[num].split()
for num1 in range(len(mylist)):
    for num2 in range(len(mylist[num1])):
        mylist[num1][num2] = int(mylist[num1][num2])

mytree = NcbiTaxonTree(1)

tempset =set()
for num1 in range(len(mylist)):
    tempset.add(mylist[num1][0])
for ele in tempset:
    mytree.add_children(NcbiTaxonTree(ele, mytree))

def add2tree(t_parent, t_child, tree):
    """
    given a tree, t_parent is the parent of t_child.
    both t_parent and t_child is int, the taxonid of NcbiTaxonTree class
    if t_parent in the tree, or sub tree of the tree, then add the t_child to the right positon, if t_child is not in tree.
    else, return information that t_parent not in the tree
    """
    if t_parent == tree.taxonid:
        if t_child not in tree.childrenids():
            tree.add_children(NcbiTaxonTree(t_child,tree))
            return True
    else:
        if tree.children() != []:
            for subtree in tree.children():
                if add2tree(t_parent, t_child,subtree):
                    break
        else:
            return str(t_parent) + "not in this tree, cannot add " + str(t_child)

def chainadd2tree(treechains, emptytree):
    """
    add2tree is too slow for large file.
    actually, the tree chain is more important. for add2tree, each time the whole tree is visited. 
    with treechain, only those branches needed were visited.
    treechain is like "[131567,2759,33154,33208,6072]"
    first adding all leaf to the tree, allowing dupliate taxonid
    then remove the duplicated ones
    treechains is a list of list. the final leaf is at the end of the tree.
    """
    for num1 in range(len(treechains)):
        temptree = emptytree
        for num2 in range(len(treechains[num1])):
            temptree.add_children(NcbiTaxonTree(treechains[num1][num2], temptree))
            temptree = temptree.children()[-1]
    return emptytree

def combinebranch(mytree):
    """
    for tree generated in chianadd2tree, there are a lot of branches with the same nodes
    return a new tree after comine those same nodes
    """
    sublist = mytree.children()
    if sublist == []:
        return mytree
    elif len(sublist) == 1:
        return mytree
    else:
        tempset = set(mytree.childrenids())
        mytree._children =[]
        for ele in tempset:
            mytree.add_children(NcbiTaxonTree(ele, mytree))
        for subtree1 in sublist:
            for subtree2 in mytree.children():
                if subtree1.taxonid == subtree2.taxonid and subtree1.children() != []:
                    subtree2.add_children(subtree1.children()[0])
        for subtree2 in mytree.children():
            combinebranch(subtree2)
        return mytree

mytree= NcbiTaxonTree(1)
mytree = chainadd2tree(mylist,mytree)
combinebranch(chainadd2tree(mylist,mytree))

for num1 in range(len(mylist)):
    add2tree(1,mylist[num1][0], mytree)
    print(num1)
    for num2 in range(len(mylist[num1])-1):
        add2tree(mylist[num1][num2], mylist[num1][num2 + 1], mytree)
            
            

