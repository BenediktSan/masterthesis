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

import header_constants as hc



if os.path.exists("build") == False:
    os.mkdir("build")

if os.path.exists("build/plots") == False:
    os.mkdir("build/plots")

if os.path.exists("build/plots/averaged_correlation") == False:
    os.mkdir("build/plots/averaged_correlation")

#include constants to rescale time

#constants from book 
gamma = 251.662 * 10**6 # rad * µs/T   40.077;//40.069244;//40.077 ;// 1/(T µs)
lattice_const = 2.7315 *10**-10# 5.451;// 2.724;//5.451; //in Angstrom

hbar =const.hbar
mu_0 = const.mu_0

coupling = mu_0 * gamma**2 * hbar**2 / ( 4 * np.pi *lattice_const**3)
factor = hbar**2  * 10**30 
#print(f"Coupling {coupling}\nfactor {factor}\nµ_0 {mu_0}\nhbar {hbar}")

##crawl through plot folders in master thesis folder to replace plots if ithey were handpicked to be in the thesis
# only does this if at runtime "renew_all_plots" is given as argument
#otherwise returns empyt list

def find_pdfs( target_filename: str, root_folder: str = "../masterthesis/images/plots"):
    """
    Recursively searches for PDF files with the given name in all subfolders.

    Parameters:
        root_folder (str): The folder to start searching from.
        target_filename (str): The exact name of the PDF file to search for (e.g., "document.pdf").

    Returns:
        list[str]: A list of full paths to matching PDF files.
    """
    matching_paths = []

    # Return if no arguments are given
    if len(sys.argv) < 2:
        return matching_paths

    if __name__ == "__main__":
        if sys.argv[1] == "renew_all_plots":
            print("\n\tRenewing all plots in master thesis folder")
        else:
            return matching_paths

    # Walk through all directories and subdirectories
    for dirpath, _, filenames in os.walk(root_folder):
        for file in filenames:
            if target_filename.lower() in file.lower() and file.lower().endswith(".pdf"):
                full_path = os.path.join(dirpath, file)
                matching_paths.append(full_path)

    return matching_paths



def find_data_for_comparison(filename):

    B_field_name = "[" + filename.split("[")[1].split("]")[0] + "]"

    #find starkov data

    starkov_data = h5py.File("../theoretical_data/caf2_data_starkov.hdf5", "r")

    starkov_Cx = starkov_data[B_field_name]["classical"]["Cx"]
    starkov_t = starkov_data[B_field_name]["classical"]["ts"]

    starkov_Cx_exp = starkov_data[B_field_name]["experiment"]["Cx"]
    starkov_t_exp = starkov_data[B_field_name]["experiment"]["ts"]


    #find Gräßers data

    if(B_field_name == "[001]"):
        new_B_field_name = "100"
    elif(B_field_name == "[011]"):
        new_B_field_name = "110"
    elif(B_field_name == "[111]"):
        new_B_field_name = "111"


    for file in os.listdir("../theoretical_data"):
        if ".hdf5" in file:
            if B_field_name in file and "Timo" in file and "FID" in file:
                timo_data = h5py.File("../theoretical_data/"+file, "r")

                timo_t = timo_data['CaF-' + new_B_field_name + ' FID from nl-spinDMFT'][0]
                timo_Cx = timo_data['CaF-' + new_B_field_name + ' FID from nl-spinDMFT'][5]
                break

    return starkov_Cx, starkov_t, timo_Cx, timo_t


def fit(t, C, omega, psi, gamma):
    return C * np.exp(-gamma*t) * np.cos(omega*t + psi)

def calc_fit_params(name,t, C):
    #initial guess
    C_0 = 0.6
    omega_0 = 160 #MHz
    psi_0 = -0.5
    gamma_0 = 0.052 #1/ms

    popt, pcov = curve_fit(fit, t, C, p0=[C_0, omega_0, psi_0, gamma_0])
    pcov = np.sqrt(np.diag(pcov))

    print(f"\n\tFit parameters for {name}:\n\tC = {popt[0]} {pcov[0]}\n\tomega = {popt[1] } {pcov[1]} rad/µs\n\tpsi = {popt[2]} {pcov[2]} rad\n\tgamma = {popt[3] *1000} {pcov[3] *1000} kHz\n")

    return popt, pcov




def plot_averaged_correlation(filename,FID_exp, t_exp):


    
    C, err, t= np.genfromtxt("build/daten/"+filename,skip_header= 1, unpack = True)

    #get theoretical data
    C_stark, t_stark, C_timo, t_timo = find_data_for_comparison(filename)


    #very arbtitary stuff happening here
    #factor 1/2 means i miss a factor 2 in my coupling???????????
    t_tilde = t /(hbar )*10**-30 * 10**6 
    t_tilde_2 = t_tilde /(hbar * 10**34 ) 

    #And now we also do a retrospective correction to the lattice constant, because I used a wrong one in th e simulation
    new_a = 2.72325
    my_a = 2.7315
    factor_a = (new_a/my_a)**(-3)
    print( f"factor a = {factor_a}")
    t_tilde = t_tilde / factor_a

    print(f"\t\tstepsize for t_tilde = {t_tilde[1]}µs")

    if("_old" in filename):
        t_tilde /= 2

    #Normlalization to 0
    print(f"\t\tFID(0) before normalization = {C[0]}")
    C = C / C[0]

    #mask to have starkovs data end when the experimental data ends
    mask_starkov = t_stark < t_exp[-1]
    C_stark = C_stark[mask_starkov]
    t_stark = t_stark[mask_starkov]

    #C_stark= np.genfromtxt("../theoretical_data/Cx.txt", unpack = True)
    #t_stark= np.genfromtxt("../theoretical_data/ts.txt", unpack = True)
#
    #C_timo= np.genfromtxt("../theoretical_data/FID_timo.txt", unpack = True)
    #t_timo= np.genfromtxt("../theoretical_data/t_timo.txt", unpack = True)


    ##little experiment to check if chang the lattice constant does change something
    #starkov_a = 2.72
    #my_a = 2.7315
    #factor_a = (starkov_a/my_a)**(-3)
    #print( f"factor a = {factor_a}")
    #
    #t_tilde = t_tilde / factor_a


    params, cov = calc_fit_params(filename,t_tilde, C)

    exp_params, exp_cov = calc_fit_params("Experimental data", t_exp, FID_exp)


    size_label = hc.size_label

    plt.figure()
    plt.rc('axes', labelsize=size_label)
    plt.axhline(y = 0, xmin = 0, xmax = 1,linestyle ="dotted", color = "black")
    plt.plot([], [], color = hc.FID_color,label=r"Semi-classical FID")
    #plt.plot(t_tilde , fit(t_tilde, *params), color = hc.fit_color, linestyle ="dashed", label = r"Fit function")
    #plt.plot(t_exp , fit(t_exp, *exp_params), color = hc.exp_color, linestyle ="dashed", label = r"Experimental fit")
    #plt.plot(t_tilde , np.exp(-t_tilde * params[2]), color = hc.fit_color, linestyle ="dashed", label = r"Fit: $\exp(-\gamma t) \cos(\omega t + \varphi)$")
    #converting into µs and adding all factors neglected in simulation
    plt.plot(t_stark[:-1],C_stark[:-1], color = hc.starkov_color_FID, label=r"Starkov et al. ")
    plt.plot(t_timo,C_timo,color = hc.timo_color_FID, label=r"Gräßer et al. ")
    plt.plot(t_tilde , C, color = hc.FID_color)
    #plt.fill_between(t_tilde, C + err, C - err, color = hc.FID_color, alpha = 0.3, label ="Estimated averaging error")
    plt.plot(t_exp, FID_exp, color = hc.exp_color, label=r"Experimental measurement")
    plt.xlabel( hc.time_label )
    plt.ylabel( hc.FID_label )
    plt.legend()
    plt.tight_layout()
    pdfname = filename[:-4]
    plt.savefig("build/plots/averaged_correlation/" + pdfname + ".pdf")
    #check if plot is in master thesis folder and replace it if it is
    matching_paths = find_pdfs( pdfname)
    for path in matching_paths:
        print(f"\n\tReplacing plot in master thesis folder: {path}")
        plt.savefig(path )
    plt.close()

    return

def plot_log_correlation(filename, FID_exp, t_exp):


    C, err, t= np.genfromtxt("build/daten/"+filename,skip_header= 1, unpack = True)

    #get theoretical data
    C_stark, t_stark, C_timo, t_timo = find_data_for_comparison(filename)

    #very arbtitary stuff happening here
    #factor 1/2 means i miss a factor 2 in my coupling???????????
    t_tilde = t /(hbar )*10**-30 * 10**6 
    t_tilde_2 = t_tilde /(hbar * 10**34 ) 
    
    #And now we also do a retrospective correction to the lattice constant, because I used a wrong one in th e simulation
    new_a = 2.72325
    my_a = 2.7315
    factor_a = (new_a/my_a)**(-3)
    #print( f"factor a = {factor_a}")
    t_tilde = t_tilde / factor_a
    #print(f"\t\tmin von exp data t = {t_exp[FID_exp.argmin()]}\n\t\tmin von t_tilde t = {t_tilde[C.argmin()]}\n\t\tmin von t_tilde_2 t = {t_tilde_2[C.argmin()]}")

    if("_old" in filename):
        t_tilde /= 2

    C = C / C[0]

    #mask to have starkovs data end when the experimental data ends
    mask_starkov = t_stark < t_exp[-1]
    C_stark = C_stark[mask_starkov]
    t_stark = t_stark[mask_starkov]

    #C_stark= np.genfromtxt("../theoretical_data/Cx.txt", unpack = True)
    #t_stark= np.genfromtxt("../theoretical_data/ts.txt", unpack = True)
#
    #C_timo= np.genfromtxt("../theoretical_data/FID_timo.txt", unpack = True)
    #t_timo= np.genfromtxt("../theoretical_data/t_timo.txt", unpack = True)


    ##little experiment to check if chang the lattice constant does change something
    #starkov_a = 2.72
    #my_a = 2.7315
    #factor_a = (starkov_a/my_a)**(-3)
    #print( f"factor a = {factor_a}")
    #
    #t_tilde = t_tilde / factor_a

    size_label = hc.size_label

    plt.figure()
    plt.rc('axes', labelsize=size_label)
    plt.axhline(y = 0, xmin = 0, xmax = 1,linestyle ="dotted", color = "black")
    #converting into µs and adding all factors neglected in simulation
    plt.plot([], [], color = hc.FID_color,label=r"Semi-classical FID")
    plt.plot(t_stark[:-1],abs(C_stark[:-1]), color = hc.starkov_color_FID, label=r"Starkov et al. ")
    plt.plot(t_timo,abs(C_timo),color = hc.timo_color_FID, label=r"Gräßers et al. ")
    plt.plot(t_exp, abs(FID_exp), color = hc.exp_color, label=r"Experimental data")
    plt.plot(t_tilde , abs(C), color = hc.FID_color)
    #plt.fill_between(t_tilde, C + err, C - err, color = hc.FID_color, alpha = 0.3, label ="Estimated averaging error")
    plt.yscale("log")
    plt.xlabel( hc.time_label )
    plt.ylabel( r'$|$ ' + hc.FID_label + r'$|$' )
    plt.ylim([10**(-5),None])    
    plt.tight_layout()
    plt.legend()
    pdfname = filename[:-4] + "_log"
    plt.savefig("build/plots/averaged_correlation/" + pdfname + ".pdf")
    #check if plot is in master thesis folder and replace it if it is
    matching_paths = find_pdfs( pdfname)
    for path in matching_paths:
        print(f"\n\tReplacing plot in master thesis folder: {path}")
        plt.savefig(path )
    plt.close()


#load experimental data

exp_100 = np.array([0])
t_100 = np.array([0])

exp_110 = np.array([0])
t_110 = np.array([0])

exp_111 = np.array([0])
t_111 = np.array([0])




for file in os.listdir("../experimental_data"):
    if file.endswith(".txt"):
        #print(file)
        if "[100]" in file:
            #print(f"\n\tFile {file} opened\n")
            t_100, exp_100  = np.genfromtxt("../experimental_data/"+file, unpack = True)
            continue
        if "[110]" in file:
            #print(f"\n\tFile {file} opened\n")
            t_110, exp_110  = np.genfromtxt("../experimental_data/"+file, unpack = True)
            continue
        if "[111]" in file:
            #print(f"\n\tFile {file} opened\n")
            t_111, exp_111  = np.genfromtxt("../experimental_data/"+file, unpack = True)
            continue


#getting rid of factor 10 from data
exp_100[49::] = 1/10 * exp_100[49::]
exp_110[95::] = 1/10 * exp_110[95::]
exp_111[58::] = 1/10 * exp_111[58::]





file_amount = 0

for file in os.listdir("build/daten"):
    if file.endswith(".txt"):
        #print(file)
        if(file.endswith("averaged_Correlation.txt") == True and "[001]" in file and not "pair" in file):
            file_amount +=1        
            print(f"\n\tFile {file} opened\n")
            plot_averaged_correlation(file, exp_100, t_100)
            print("\n\tplotting log")
            plot_log_correlation(file, exp_100, t_100)
            continue
        if(file.endswith("averaged_Correlation.txt") == True and "[011]" in file and not "pair" in file):
            file_amount +=1        
            print(f"\n\tFile {file} opened\n")
            plot_averaged_correlation(file, exp_110, t_110)
            print("\n\tplotting log")
            plot_log_correlation(file, exp_110, t_110)
            continue
        if(file.endswith("averaged_Correlation.txt") == True and "[111]" in file and not "pair" in file):
            file_amount +=1        
            print(f"\n\tFile {file} opened\n")
            plot_averaged_correlation(file, exp_111, t_111)
            print("\n\tplotting log")
            plot_log_correlation(file, exp_111, t_111)
            continue

        

        
print("\nDateianzahl used for averaged system plots: ", file_amount)

