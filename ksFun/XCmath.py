# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 16:11:15 2016

@author: k
"""
import numpy as np
def gaussian(x, ampl, center, sigma):
    '''
    computes the Gaussian function
     Parameters
    ----------
    x : number
        Point to evaluate the Gaussian for.
    a : number
        Amplitude.
    b : number
        Center.
    c : number
        Width.

    Returns
    -------
    float
        Value of the specified Gaussian at *x*
    '''
    return ampl * np.exp(-(x - float(center)) ** 2 / (2.0 * sigma ** 2 ))

def gaussian_fit(x, y, plot = False):
    '''
    Performs a Gaussian fitting of the specified data.
    return amplitude, center and sigma of gaussian function
    Parameters
    ----------
    x : ndarray
        Data on the x axis.
    y : ndarray
        Data on the y axis.
    Returns
    ----------
    returns the parameters of the Gaussian that fits the specified data.
    '''
    from scipy import optimize
    ampl = np.max(y)
    center = x[0]
    sigma = (x[1] - x[0]) *5
    initial = [ampl, center, sigma]
    params, pcov = optimize.curve_fit(gaussian, x, y, initial)
    if plot:
        import matplotlib.pyplot as plt
        y2 = [gaussian(i, *params) for i in x]
        plt.plot(x,y,x,y2)
    return params

def mappingScore_geneComparison(n,a,b, ab,comparisonNum):
    """
    for hypergenometic testing in within species stage/tissue/cell comparison.
    check paper Comparison of D. melanogaster and C. elegans developmental stages, tissues, and cells by modENCODE RNA-seq data
    method part
    """
    from decimal import Decimal
    import math
    n = Decimal(n)
    a = Decimal(a)
    b = Decimal(b)
    ab = Decimal(ab)
    comparisonNum = Decimal(comparisonNum)
#    def comb(total, choose):
#        return Decimal(math.factorial(total)) / (Decimal(math.factorial(choose)) * Decimal(math.factorial(total - choose)))
#    def pvalue(i):
#        return comb(n,i) * comb(n-i,a-i) * comb(n-a, b-i) / (comb(n,a) * (comb(n,b)))
    
#    def f(value):
#        return Decimal(math.factorial(value))
#    def pvalue2(i):
#        return f(a) / f(a-i) * f(b) / f(b-i) /f(i) * f(n-a) / f(n) * f(n-b) / f(n+i-a-b)

    def fl(value):
        value = int(value)
        return Decimal(sum(math.log(i) for i in range(1,value+1)))
    def pvalue3(i):
        return Decimal(math.e) ** (fl(a) - fl(a-i) + fl(b) - fl(b-i) -fl(i) + fl(n-a) - fl(n) + fl(n-b) - fl(n+i-a-b))
    
    p = Decimal()
    for i in range(int(ab), int(min(a,b))+1):
        p += pvalue3(i)
    return - float(p.log10() + comparisonNum.log10())
    
    
    
