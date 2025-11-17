import numpy as np
import matplotlib.pyplot as plt
#import uncertainties as unc
#import uncertainties.unumpy as unp
#from uncertainties import ufloat
from scipy.optimize import curve_fit
import scipy.constants as const
#import sympy
import os
import re
import sys



if os.path.exists("build") == False:
    os.mkdir("build")

if os.path.exists("build/plots2D") == False:
    os.mkdir("build/plots2D")

if os.path.exists("build/plots2D/variance") == False:
    os.mkdir("build/plots2D/variance")

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


def extract_number(s):
    match = re.search(r'\d+', s)
    return int(match.group()) if match else None






def plot_variances_for_certain_times(variance_filename, FID_filename):#,FID_exp, t_exp):

    print(f"\nstarted plotting for \n{variance_filename}\n\n")


    old_name = ""
    if "_old" in variance_filename:
        old_name = "old"
    

    #read data
    variance, exp_val, exp_val_squared, t_var = np.genfromtxt("build/daten_2D/"+variance_filename,skip_header= 1, unpack = True)
    C, err, t= np.genfromtxt("build/daten_2D/"+FID_filename,skip_header= 1, unpack = True)

    #normalize in case its not normalized

    C_norm = C / C[0]

    #recaclulate variance
    exp_val_norm = exp_val / exp_val[0]
    exp_val_squared_norm = exp_val_squared / exp_val[0]**2
    variance_norm = exp_val_squared_norm - exp_val_norm**2
    #print(variance_norm)

 
    size_label = 15

    #print(f"\t\tsize of t_tilde {t_tilde.size} and size of pair_corrs {pair_correlations[0].size} {pair_correlations[1].size}")


    plt.figure()
    plt.rc('axes', labelsize=size_label)
    #converting into µs and adding all factors neglected in simulation
    plt.plot(t , C_norm, color = "green",label=r"$C_x$")
    for i in range(0, len(t_var)):
        plt.plot(t_var[i] , exp_val_norm[i], color = "green")
        plt.errorbar(t_var[i], exp_val_norm[i], yerr = np.sqrt(variance_norm[i]), fmt ="^", color = plt.gca().lines[-1].get_color(), capsize=3, label=r"Variance for "+ str(round(t_var[i])) + r" /J")
    plt.xlabel( r'$t $ / $J$')
    plt.ylabel( r'Var$[FID]$ ' )
    plt.legend()
    pdfname = variance_filename[:-4]
    plt.savefig("build/plots2D/variance/" + pdfname + ".pdf")
    #check if plot is in master thesis folder and replace it if it is
    matching_paths = find_pdfs( pdfname)
    for path in matching_paths:
        print(f"\n\tReplacing plot in master thesis folder: {path}")
        plt.savefig(path )    
    plt.close()

    return

def plot_ALL_variances_for_certain_times(grouped_filenames_variance, FID_filenames, B_field_name):

    old_name = []

    for file in grouped_filenames_variance:
        if "_old" in file:
            old_name.append("old")
        else:
            old_name.append("")

    old_name = np.asarray(old_name)




    if grouped_filenames_variance.size == 1:
        return
    print(f"\nstarted plotting for \n{grouped_filenames_variance}\n\n")

    t_values =[]
    variance =[]
    exp_val  =[]
    exp_val_squared = []

    for file in grouped_filenames_variance:
        #read data
        variance_temp, exp_val_temp, exp_val_squared_temp, t_values_temp = np.genfromtxt("build/daten_2D/"+file,skip_header= 1, unpack = True)
        t_values.append(t_values_temp)
        variance.append(variance_temp)
        exp_val.append(exp_val_temp)
        exp_val_squared.append(exp_val_squared_temp)
    

    t_values = np.asarray(t_values)
    variance = np.asarray(variance)
    exp_val = np.asarray(exp_val)
    exp_val_squared = np.asarray(exp_val_squared)


    #red FID data
    t = []
    C = []
    err = []


    for file in FID_filenames:
        C_temp, err_temp, t_temp = np.genfromtxt("build/daten_2D/"+file,skip_header= 1, unpack = True)
        t.append(t_temp)
        C.append(C_temp)
        err.append(err_temp)

    t = np.asarray(t)
    C = np.asarray(C)
    err = np.asarray(err)

    #normalize in case its not normalized

    for i in range(0, len(C)):
        norm = exp_val[i][0]

        C[i] = C[i] / C[i][0]
        exp_val[i] = exp_val[i] / norm
        exp_val_squared[i] = exp_val_squared[i] / norm**2

    #recaclulate variance

    new_variance_norm = []
    for i in range(0, len(variance)):
        new_variance_norm.append(exp_val_squared[i] - exp_val[i]**2)

    new_variance_norm = np.asarray(new_variance_norm)
    #print(new_variance_norm)


  
    time = []
    time_var = []

    for i in range(0, len(new_variance_norm)):
        
        time.append(t[i])
        time_var.append(t_values[i])

    time = np.asarray(time)
    time_var = np.asarray(time_var)



    size_label = 15

    plt.figure()
    plt.rc('axes', labelsize=size_label)
    #converting into µs and adding all factors neglected in simulation
    for i in range(0, len(grouped_filenames)):
        special_case = ""
        if "no_norm" in grouped_filenames[i]:
            special_case = " no norm"
        plt.plot(time_var[i], exp_val[i], "^", label= "Variance for "+ str(extract_number(grouped_filenames[i][20:])) + r"$^2$ Particles" + special_case)
        plt.errorbar(time_var[i], exp_val[i], yerr = np.sqrt(new_variance_norm[i]), fmt ="^", color = plt.gca().lines[-1].get_color(), capsize=3)
        plt.plot(time[i] , C[i], color = plt.gca().lines[-1].get_color())

    plt.xlabel( r'$t $ / $J$')
    plt.ylabel( r'$Var[FID]$ ' )
    plt.legend()
    pdfname = "2D_variance_comparison_" + B_field_name + "_" + str(extract_number(grouped_filenames[0][3:]))
    plt.savefig("build/plots2D/variance/" + pdfname + ".pdf")
    #check if plot is in master thesis folder and replace it if it is
    matching_paths = find_pdfs( pdfname)
    for path in matching_paths:
        print(f"\n\tReplacing plot in master thesis folder: {path}")
        plt.savefig(path )    
    plt.close()



all_filenames = np.array([])

for file in os.listdir("build/daten_2D"):
    all_filenames = np.append(all_filenames,file)

#get rid of pair correlation files



#mask = "pair" not in all_filenames
#all_filenames = all_filenames[mask]
mask = [False] * len(all_filenames)
for i in range(0, len(all_filenames)):
    mask[i] = "variance"  in all_filenames[i] and "2D" in all_filenames[i]
all_variance_filenames = all_filenames[mask]

allready_checked_sys_amount = np.array([])
allready_checked_B_field = np.array([])



for file in all_variance_filenames:

    #sort for B_fields
    B_field_name = "[" + file.split("[")[1].split("]")[0] + "]"

    
    system_amount = extract_number(file[3:])

    continue_flag = True

    for i in range(0, allready_checked_sys_amount.size):
        if system_amount == allready_checked_sys_amount[i] and B_field_name == allready_checked_B_field[i]:
            #print(f"Already checked {system_amount} and {B_field_name}")
            continue_flag = False

    if continue_flag == False:
        continue

    mask = [False] * len(all_variance_filenames)
    for i in range(0, len(all_variance_filenames)):
        if(system_amount == extract_number(all_variance_filenames[i][3:]) and B_field_name in all_variance_filenames[i]):
            mask[i] = True
    grouped_filenames = all_variance_filenames[mask]

    #ffind files for FID signal
    FID_filenames = np.array([])
    
    for val_file in grouped_filenames:
        val_file_name = str(val_file).split("_variance")[0]
        for i in range(0, len(all_filenames)):
            if val_file_name in all_filenames[i] and "variance" not in all_filenames[i] and "pair" not in all_filenames[i] and "interim" not in all_filenames[i] and "copy" not in all_filenames[i]:
                FID_filenames =  np.append(FID_filenames, all_filenames[i])
                


    allready_checked_sys_amount = np.append(allready_checked_sys_amount, system_amount)
    allready_checked_B_field = np.append(allready_checked_B_field, B_field_name)

    if grouped_filenames.size != FID_filenames.size:
        print(f"\n\tWARNING: Different amount of grouped and FID files for system amount {system_amount}.\n")
        print(f"\t\tThis case is skipped.\n")
        continue

    if grouped_filenames.size == 0:
        continue

    
    plot_ALL_variances_for_certain_times(grouped_filenames, FID_filenames, B_field_name)#, exp_111, t_111)                    plot_interim_difference(file, interim_results_file)
    for i in range(0, grouped_filenames.size):
        plot_variances_for_certain_times(grouped_filenames[i], FID_filenames[i])
