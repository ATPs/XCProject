# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 19:37:26 2017

@author: ATPs
"""

def basicPractice20170902():
    from ete3 import Tree
    
    def combineTree(t1,t2,align_t1, align_t2):
        '''
        combine Tree
        '''
        if t1.is_leaf() and t2.is_leaf():
            if t1.name in align_t1 and t2.name in align_t2:
                if align_t1.index(t1.name) == align_t2.index(t2.name):
                    return Tree(name = ';'.join([t1.name,t2.name]))
                else:
                    print("Not the same position")
            elif t1.name not in align_t1 and t2.name not in align_t2:
                t = Tree()
                t.add_child(t1)
                t.add_child(t2)
                return t
            else:
                print('something wrong?')
        else:
            print('wait')
        
    t1 = Tree(name = 'A1')
    t2 = Tree(name = 'A2')
    align_t1 = []
    align_t2 = []
    t = combineTree(t1,t2,align_t1, align_t2)
    print(t)

