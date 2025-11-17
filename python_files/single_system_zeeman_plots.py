import numpy as np
import matplotlib.pyplot as plt
#import uncertainties as unc
#import uncertainties.unumpy as unp
#from uncertainties import ufloat
from scipy.optimize import curve_fit
import scipy.constants as const
#import sympy
import os


if os.path.exists("build") == False:
    os.mkdir("build")

if os.path.exists("build/plots") == False:
    os.mkdir("build/plots")






def plot_magnetization(filename):
    
    Mx, My, Mz, t= np.genfromtxt("build/daten/"+filename,skip_header= 1, unpack = True)

    size_label = 15

    plt.figure()
    plt.rc('axes', labelsize=size_label)
    plt.plot(t, Mx, label=r"M_x")
    plt.plot(t, My, label=r"M_y")
    plt.plot(t, Mz, label=r"M_z")
    plt.plot(t, np.sqrt(Mx**2 + My**2 + Mz**2),label=r"|M|")
    plt.plot(t, np.sqrt(Mx**2 + My**2),label=r"|M_{xy}|")
    plt.xlabel( r'$t$')
    plt.ylabel( r'$M$ ' )
    plt.legend()
    plt.savefig("build/plots/" + filename[:-4] +".pdf")
    plt.close()

    return

def plot_magnetization_CS(filename):
    
    Sx, Sy, Sz, t= np.genfromtxt("build/daten/"+filename,skip_header= 1, unpack = True)

    size_label = 15

    plt.figure()
    plt.rc('axes', labelsize=size_label)
    plt.plot(t, Sx, label=r"S_x")
    plt.plot(t, Sy, label=r"S_y")
    plt.plot(t, Sz, label=r"S_z")
    plt.plot(t, np.sqrt(Sx**2 + Sy**2 + Sz**2),label=r"|S|")
    plt.xlabel( r'$t$')
    plt.ylabel( r'$S$ ' )
    plt.legend()
    plt.savefig("build/plots/" + filename[:-4] +".pdf")
    plt.close()

    return

def plot_correlation(filename):
    
    G, t= np.genfromtxt("build/daten/"+filename,skip_header= 1, unpack = True)

    size_label = 15

    plt.figure()
    plt.rc('axes', labelsize=size_label)
    plt.plot(t, G, label=r"G")
    plt.xlabel( r'$t$')
    plt.ylabel( r'$G$ ' )
    plt.legend()
    plt.savefig("build/plots/" + filename[:-4] +".pdf")
    plt.close()

    return



file_amount = 0

for file in os.listdir("build/daten"):
    if file.endswith(".txt"):
        #print(file)
        if(file.endswith("Zeeman_Magnetization.txt") ==True ):
            file_amount +=1        
            print(f"\n\tFile {file} opened\n")
            plot_magnetization(file)
            continue
        if(file.endswith("Zeeman_Magnetization_CS.txt") ==True ):
            file_amount +=1        
            print(f"\n\tFile {file} opened\n")
            plot_magnetization_CS(file)
            continue
        if(file.endswith("Zeeman_Correlation.txt") ==True ):
            file_amount +=1        
            print(f"\n\tFile {file} opened\n")
            plot_correlation(file)
            continue
        
        

        
print("Dateianzahl opened for single system Zeeman plots: ", file_amount)

