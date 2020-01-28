# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 15:24:58 2019

@author: ATPs
"""

def merge(intervals):
    '''
    for a list of intervals, merge those overlap with each other
    '''
    interv=intervals.copy()
    interv.sort()
    res=[]
    while(len(interv)>0):
        if len(interv)==1:
            res.append(interv[0])
            interv.pop(0)
            continue
        if interv[0][1]>=interv[1][0]:
            tmp=[interv[0][0],max(interv[0][1],interv[1][1])]
            interv[0]=tmp
            interv.pop(1)
            continue
        res.append(interv[0])
        interv.pop(0)
    return res