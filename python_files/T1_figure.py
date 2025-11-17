import numpy as np
import matplotlib.pyplot as plt
#import uncertainties as unc
#import uncertainties.unumpy as unp
#from uncertainties import ufloat
from scipy.optimize import curve_fit
import scipy.constants as const
#import sympy
import os
import h5py 
import sys

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command


import header_constants as hc


def T1_func(t, ton, T1):
    func = np.zeros(len(t))
    ton = 2.0

    for i in range(len(t)):
        if t[i] < ton:
            func[i] = 0
        else:
            func[i] = 1 - np.exp(-(t[i] - ton) / T1)
    return func

Tmax = 14.0
t= np.linspace(0, Tmax, 1000)
T_1 = 2.0
ton = 2.0

size_label = hc.size_label


fig, ax = plt.subplots()

ax.axhline(y = 0, xmin = 0, xmax = 1,linestyle ="dotted", color = "black")
ax.plot(t , T1_func(t, ton, T_1), color = "teal")
ax.set_ylabel( r"$M^z$",fontsize=size_label)
ax.set_xlabel( "time" ,fontsize=size_label)
ax.set_yticks([1])
ax.set_yticklabels([r"$M_\text{equi}^z$"], fontsize=14)
ax.set_xticks([ton])
ax.set_xticklabels([r"$t_\text{on}$"], fontsize=14)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
#plt.legend()
fig.tight_layout()
plt.savefig("../masterthesis/images/images/T1_figure.pdf")