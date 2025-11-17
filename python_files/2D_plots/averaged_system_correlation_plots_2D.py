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

sys.path.append('/home/sander/Data/Masterarbeit/python_files/')

sys.path.append('/mnt/d/Dateien/Uni/Masterarbeit/python_files/')
import header_constants as hc



if os.path.exists("build") == False:
    os.mkdir("build")

if os.path.exists("build/plots2D") == False:
    os.mkdir("build/plots2D")

if os.path.exists("build/plots2D/averaged_correlation") == False:
    os.mkdir("build/plots2D/averaged_correlation")

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

def find_pdfs( target_filename: str, root_folder: str = "../masterthesis/images/plots_2D"):
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


    #find timos data

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



def plot_averaged_correlation(filename):

    old_name = ""
    if("_old" in filename):
        old_name = "_old"
    
    C, err, t= np.genfromtxt("build/daten_2D/"+filename,skip_header= 1, unpack = True)

    #get theoretical data
    C_stark, t_stark, C_timo, t_timo = find_data_for_comparison(filename)


    #normalization

    #C = C / C[0]

    #C_stark= np.genfromtxt("../theoretical_data/Cx.txt", unpack = True)
    #t_stark= np.genfromtxt("../theoretical_data/ts.txt", unpack = True)
#
    #C_timo= np.genfromtxt("../theoretical_data/FID_timo.txt", unpack = True)
    #t_timo= np.genfromtxt("../theoretical_data/t_timo.txt", unpack = True)


    size_label = hc.size_label

    plt.figure()
    plt.rc('axes', labelsize=size_label)
    plt.axhline(y = 0, xmin = 0, xmax = 1,linestyle ="dotted", color = "black")
    #converting into µs and adding all factors neglected in simulation
    plt.plot(t , C, color = "blue",label=r"$C_x$")
    plt.fill_between(t, C + err, C - err, color = "blue", alpha = 0.3, label ="Estimated averaging error")
    #plt.plot(t_2 , C, label=r"C3", alpha = 0.5)
    #plt.plot(t_timo,C_timo,color = "green", label=r"Timos data ")
    plt.xlabel(hc.time_label_2D)
    plt.ylabel( r'$FID$ ' )
    plt.legend()
    pdfname = filename[:-4]
    plt.savefig("build/plots2D/averaged_correlation/" + pdfname + ".pdf")
    #check if plot is in master thesis folder and replace it if it is
    matching_paths = find_pdfs( pdfname)
    for path in matching_paths:
        print(f"\n\tReplacing plot in master thesis folder: {path}")
        plt.savefig(path )    
    plt.close()

    return





file_amount = 0

for file in os.listdir("build/daten_2D"):
    if file.endswith(".txt"):
        #print(file)
        if(file.endswith("averaged_Correlation.txt") == True and "2D" in file and "[001]" in file and not "pair" in file):
            file_amount +=1        
            print(f"\n\tFile {file} opened\n")
            plot_averaged_correlation(file)
            continue
        if(file.endswith("averaged_Correlation.txt") == True and "2D" in file and "[011]" in file and not "pair" in file):
            file_amount +=1        
            print(f"\n\tFile {file} opened\n")
            plot_averaged_correlation(file)
            continue
        if(file.endswith("averaged_Correlation.txt") == True and "2D" in file  and "[111]" in file and not "pair" in file):
            file_amount +=1        
            print(f"\n\tFile {file} opened\n")
            plot_averaged_correlation(file)
            continue

        

        
print("Dateianzahl used for averaged system plots: ", file_amount)

