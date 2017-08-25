from scipy.optimize import curve_fit
def v(s,Vmax,Km):
    return Vmax *s/(Km+s)
import numpy as np
import matplotlib.pyplot as plt

s = np.array([0,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000])
v_ob =np.array( [0,0.011,0.021,0.034,0.048,0.056,0.068,0.08,0.09,0.103,0.12,0.242,0.321,0.346,0.389,0.415,0.455,0.493,0.517,0.539])

popt, pcov=curve_fit(v,s,v_ob)
v_predict = np.array([v(e,popt[0],popt[1]) for e in s])

plt.plot(s,v_ob,s,v_predict)