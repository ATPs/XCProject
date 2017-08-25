# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 13:36:43 2016

@author: k
"""
import multiprocessing
from multiprocessing import Pool

def f(x):
    return x*x

if __name__ == '__main__':
    with Pool(5) as p:
        print(p.map(f, [1, 2, 3]))


from multiprocessing import Pool
pool = Pool()
result = pool.apply_async(f,[1,2])
answer = result.get(timeout=1)