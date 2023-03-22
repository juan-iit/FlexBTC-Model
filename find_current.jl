#************************************************************************
# Packages
# using Conda
# Conda.add("scipy")
using PyCall
#************************************************************************
# Python code
py"""
#************************************************************************
import numpy as np
from scipy.optimize import fsolve
#************************************************************************
a1	= 1.5184
a2	= 1.5421E-03
a3	= 9.523E-05
a4	= 9.84E-08
r1 = 4.45153E-05
r2 = 6.88874E-09
d1 = -3.12996E-06
d2 = 4.47137E-07
s = 0.33824
t1 = -0.01539
t2 = 2.00181
t3 = 15.24178
#************************************************************************
#Functions
def U_rev(temp):
    T = temp + 273.15
    U_rev = a1 - a2 * T + a3 * T * np.log(T) + a4 * T**2
    return U_rev

def find_i_from_p(p_val,C_E,n_cell,A_cell,T,p):
    i_pval = np.zeros(len(p_val))
    for j in range(len(p_val)):
        P_j=p_val[j] * C_E * 10**6
        def f(x):
            return (U_rev(T)+(r1+d1+r2*T+d2*p)*x+s*np.log10((t1+t2/T+t3/T**2)*x+1))*x*n_cell*A_cell - P_j
        i_pval[j] = fsolve(f,[2000])
    return i_pval
"""
#************************************************************************