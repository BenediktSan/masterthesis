import numpy as np
import matplotlib.pyplot as plt
#import uncertainties as unc
#import uncertainties.unumpy as unp
#from uncertainties import ufloat
from scipy.optimize import curve_fit
import scipy.constants as const
#import sympy
import re
import os
import h5py
import sys

import header_constants as hc   


if os.path.exists("build") == False:
    os.mkdir("build")

if os.path.exists("build/plots") == False:
    os.mkdir("build/plots")


if os.path.exists("build/plots/averaged_pair_correlation") == False:
    os.mkdir("build/plots/averaged_pair_correlation")


#include constants to rescale time

#constants from book 
gamma = 251.662 * 10**6 # rad * µs/T   40.077;//40.069244;//40.077 ;// 1/(T µs)
lattice_const = 2.7315 *10**-10# 5.451;// 2.724;//5.451; //in Angstrom

hbar =const.hbar
mu_0 = const.mu_0

coupling = mu_0 * gamma**2 * hbar**2 / ( 4 * np.pi *lattice_const**3)
factor = hbar**2  * 10**30 
#print(f"Coupling {coupling}\nfactor {factor}\nµ_0 {mu_0}\nhbar {hbar}")




####Functions

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



def read_txt_file(file_path, delimiter=None):
    """
    Reads a text file dynamically without knowing the number of columns.
    Automatically detects the delimiter if not provided.
    
    :param file_path: Path to the text file
    :param delimiter: Delimiter used in the file (e.g., ',', '\t', ' ')
    :return: List of lists containing the data
    """
    data = []
    
    with open(file_path, 'r') as file:
        skip_lines = 1
        for _ in range(skip_lines):
            next(file, None)  # Skip the specified number of lines
        

        for line in file:
            line = line.strip()
            if line:
                # Detect delimiter if not provided
                if delimiter is None:
                    delimiter = '\t' if '\t' in line else ',' if ',' in line else ' '
                
                # Split the line using the detected delimiter
                columns = line.split(delimiter)
                data.append(columns)
    
    # Transpose the data to get columns instead of rows
    data = list(map(list, zip(*data)))

    return data


def extract_number(s):
    match = re.search(r'\d+', s)
    return int(match.group()) if match else None

##functions for interpolation

def interpolate(t_value, t_exp, exp_data):

    #special case if t is larger than my maximal experimental time
    if t_value > t_exp[-1]:

        #just retruning last value
        return exp_data[-1]
    #find values between whicht t_value lies
    mask1 = t_value >= t_exp
    mask2 = t_value < t_exp

    t_1 = t_exp[mask1][-1]
    t_2 = t_exp[mask2][0]

    data1 = exp_data[mask1][-1]
    data2 = exp_data[mask2][0]
    interp_value = np.interp(t_value, [t_1, t_2], [data1, data2])
    
    #print(f"Interpolation at {t_value} yields {interp_value}\nCalculated from ({t_1} {data1}) and ({t_2} {data2})")
    
    return interp_value

def find_data_for_comparison(filename):

    B_field_name = "[" + filename.split("[")[1].split("]")[0] + "]"

    #find Gräßers data
    new_B_field_name = ""

    if(B_field_name == "[001]"):
        new_B_field_name = "100"
    elif(B_field_name == "[011]"):
        new_B_field_name = "110"
    elif(B_field_name == "[111]"):
        new_B_field_name = "111"



    for file in os.listdir("../theoretical_data"):
        if ".hdf5" in file:
            if B_field_name in file and "Timo" in file and "pair" in file:
                timo_data = h5py.File("../theoretical_data/"+file, "r")

                timo_t = timo_data['CaF-' + new_B_field_name +' p.-c. contr. G^xx from nl-spinDMFT'][0]
                timo_auto_correlation = timo_data['CaF-' + new_B_field_name +' p.-c. contr. G^xx from nl-spinDMFT'][1]
                timo_pair_correlations = timo_data['CaF-' + new_B_field_name +' p.-c. contr. G^xx from nl-spinDMFT'][2:-1]
                break

    return timo_t, timo_auto_correlation, timo_pair_correlations



def plot_pair_comparison(pair_filename):

    old_name = []
    if "_old" in pair_filename:
        old_name = "_old"
    else:
        old_name = ""


    dimension_name = []

    if("2D" in pair_filename):
        dimension_name = "2D_"
    else:  
        dimension_name = "3D_"

    print(f"\nstarted plotting for \n{pair_filename}\n\n")


    B_field_name =  "[" + pair_filename.split("[")[1].split("]")[0] + "]"

    #find Gräßers data


    system_size = extract_number(pair_filename[10:30])
    system_size_name = str(system_size) + r"P"
    print(f"\t\tSystem Size: {system_size} Particles")

    #get pair correlations

    ##array for pair correlations
    #set to 0 for all correlations
    #this cutoff defines the amunt of correlations per subplot
    #this excludes the auto correlation plot
    correlation_cut_off = 3
    sum_of_all_correlations =  np.array([])
    t_pair = np.array([])

    pair_correlations = np.array([])

    data = read_txt_file("build/daten/"+pair_filename)  

    if(len(data) == 0):
        return  
    
    #convert to np.array
    pair_correlations = np.array([data[0]], dtype=np.float32)
    for i in range(1, len(data)):
        pair_correlations = np.append(pair_correlations,np.array([data[i]],dtype=np.float32), axis = 0)


    #extract t
    t_pair = pair_correlations[-1]
    pair_correlations = np.delete(pair_correlations,-1,0)

#calculate sum of pair correlations
    sum_of_all_correlations = pair_correlations[0]
    for i in range(1, len(pair_correlations)):
        sum_of_all_correlations = sum_of_all_correlations + pair_correlations[i]
    print(f"\t\tsum of all correalations = {sum_of_all_correlations[0]}\n")

    print(f"\t\t{len(pair_correlations)} Different Pair Correlations found")
    #extract auto correlation

    #reorder some correlations for B011
    #if B_field_name == "[011]":
    #    print(f"\t\tReordering pair correlations for B field {B_field_name}")
    #    temp1 = pair_correlations[6].copy()
    #    temp2 = pair_correlations[7].copy()
    #    temp3 = pair_correlations[8].copy()
    #    temp4 = pair_correlations[9].copy()
#
    #    pair_correlations[6] = temp2
    #    pair_correlations[9] = temp1
    #    pair_correlations[7] = temp3
    #    pair_correlations[8] = temp4

    if B_field_name == "[011]":
        print(f"\t\tReordering pair correlations for B field {B_field_name}")

        pair_correlations[6] += pair_correlations[7]
        pair_correlations = np.delete( pair_correlations,7,0)



    auto_correlation = pair_correlations[0]
    pair_correlations = np.delete(pair_correlations,0,0)


    if (len(pair_correlations) < correlation_cut_off * 3):
        print(f"\t\tWARNING: Not enough pair correlations found for {pair_filename} to plot all subplots")
        for i in range(len(pair_correlations), correlation_cut_off * 3):
            pair_correlations = np.append(pair_correlations, np.zeros((1, len(t_pair))), axis=0)
        print(f"\t\t\tAdding {(3 * correlation_cut_off) - len(pair_correlations)} empty pair correlation  to fill up subplots")


    #very arbtitary stuff happening here
    #factor 1/2 means i miss a factor 2 in my coupling???????????
    t_tilde = t_pair /(hbar )*10**-30 * 10**6 
        
    #And now we also do a retrospective correction to the lattice constant, because I used a wrong one in th e simulation
    new_a = 2.72325
    my_a = 2.7315
    factor_a = (new_a/my_a)**(-3)
    #print( f"factor a = {factor_a}")
    t_tilde = t_tilde / factor_a
    
    #print(f"\t\tmin von exp data t = {t_exp[FID_exp.argmin()]}\n\t\tmin von t_tilde t = {t_tilde[C.argmin()]}\n\t\tmin von t_tilde_2 t = {t_tilde_2[C.argmin()]}")
    #print(f"\n\t\tstepsize for t_tilde = {t_tilde[1]}")

    if "old" in pair_filename:
        t_tilde = t_tilde / 2

    ##normalizing
    
    norm =  sum_of_all_correlations[0]

    #norm = 5 ** 3 * 8
    sum_of_all_correlations = sum_of_all_correlations / norm
    auto_correlation = auto_correlation / norm
    pair_correlations = pair_correlations / norm

    ##extract Gräßers data from the hdf5 file

    timo_t, timo_auto, timo_pair = find_data_for_comparison(pair_filename)

    print(f"\t\tTimo data has {len(timo_pair)} pair correlations and {len(timo_t)} time steps")
    size_label = hc.size_label

    #print(f"\t\tsize of t_tilde {t_tilde.size} and size of pair_corrs {pair_correlations[0].size} {pair_correlations[1].size}")
    pastel1 = plt.get_cmap('Pastel1')
    pastel2 = plt.get_cmap('Pastel2')
    pastel1.colors += pastel2.colors


    fig, axs = plt.subplots(2, 2, sharex=True, figsize=(10, 8))
    #fig.suptitle(r"Pair Correlation Comparison " + str(extract_number(pair_filename[:8])) + r" i.c. and " +str(extract_number(pair_filename[10:])) +  r"$^3$ Particles")
    axs[0, 0].plot(t_tilde, auto_correlation, label =r"Semi-classical $m_c G_0^{xx}$")
    axs[0, 0].plot(timo_t,timo_auto, color = axs[0, 0].lines[-1].get_color(), linestyle = "--", label = r"nl-spinDMFT $G_0^{xx}$")
    axs[0,0].legend()
    axs[0, 0].set_title('Auto Correlation')
    for i in range(0, correlation_cut_off):
        #factor 1 due to extracted auto correlation
        c = str(i + 1)
        axs[0, 1].plot(t_tilde, pair_correlations[i], label = r"Semi-classical $m_c G_" + f"{c}" + r"^{xx}$")
        axs[0, 1].plot(timo_t, timo_pair[i], color = axs[0, 1].lines[-1].get_color(), linestyle = "--", label = r"nl-spinDMFT $m_" + f"{c}" + r" G_" + f"{c}" + r"^{xx}$")
    axs[0, 1].legend()
    #axs[0, 1].set_title('B')
    for i in range(correlation_cut_off, 2*correlation_cut_off):
        #factor 1 due to extracted auto correlation
        c = str(i + 1)
        axs[1, 0].plot(t_tilde,pair_correlations[i], label = r"Semi-classical $m_c G_{" + f"{c}" + r"}^{xx}$")
        axs[1, 0].plot(timo_t,timo_pair[i], color = axs[1, 0].lines[-1].get_color(), linestyle = "--", label = r"nl-spinDMFT $m_{" + f"{c}"  + r"} G_{" + f"{c}" + r"}^{xx}$")
    axs[1,0].legend()
    for i in range(2*correlation_cut_off, 3*correlation_cut_off):
        #factor 1 due to extracted auto correlation
        c = str(i + 1)
        axs[1, 1].plot(t_tilde, pair_correlations[i], label = r"Semi-classical $m_c G_{" + f"{c}" + r"}^{xx}$")
        axs[1, 1].plot(timo_t,timo_pair[i], color = axs[1, 1].lines[-1].get_color(), linestyle = "--", label = r"nl-spinDMFT $m_{" + f"{c}" + r"} G_{" + f"{c}" + r"}^{xx}$")
    axs[1, 1].legend()
    #axs[1, 1].set_title('C')

    #for ax in axs.flat:
    #    ax.set(xlabel=r'$t$ / $µs$', ylabel=r'$FID$ ')

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    #for ax in axs.flat:
    #    ax.label_outer()
    axs[1, 0].set_xlabel(hc.time_label,fontsize=size_label)
    axs[1, 1].set_xlabel(hc.time_label,fontsize=size_label)
    axs[0, 0].set_ylabel(hc.pair_correlation_label,fontsize=size_label)
    axs[1, 0].set_ylabel(hc.pair_correlation_label,fontsize=size_label)
    fig.tight_layout()

    pdfname = dimension_name + system_size_name + "_" + str(extract_number(pair_filename[:8])) +"sys_" + B_field_name
    plt.savefig("build/plots/averaged_pair_correlation/" + pdfname + ".pdf")
    #check if plot is in master thesis folder and replace it if it is
    matching_paths = find_pdfs( pdfname)
    for path in matching_paths:
        print(f"\n\tReplacing plot in master thesis folder: {path}")
        plt.savefig(path )
    plt.close()
    return




all_filenames = np.array([])

for file in os.listdir("build/daten"):
    all_filenames = np.append(all_filenames,file)



mask = [False] * len(all_filenames)
for i in range(0, len(all_filenames)):
    mask[i] = "pair" in all_filenames[i] and "interim" not in all_filenames[i] and "variance" not in all_filenames[i]
all_filenames = all_filenames[mask]

# just so its not empty and i can append
allready_checked_sys_amount = np.array([])
allready_checked_B_field = np.array([])


for file in all_filenames:
    
    system_amount = extract_number(file[:5])
    #sort for B_fields
    B_field_name = "[" + file.split("[")[1].split("]")[0] + "]"

    for i in range(0, allready_checked_sys_amount.size):
        if system_amount == allready_checked_sys_amount[i] and B_field_name == allready_checked_B_field[i]:
            #print(f"Already checked {system_amount} and {B_field_name}")
            break

    mask = [False] * len(all_filenames)
    for i in range(0, len(all_filenames)):
        if(system_amount == extract_number(all_filenames[i][:8]) and B_field_name in all_filenames[i]):
            mask[i] = True
    grouped_filenames = all_filenames[mask]

    #print(f"system amount {system_amount} mask{mask} groups{grouped_filenames}")

    allready_checked_sys_amount = np.append(allready_checked_sys_amount, system_amount)
    allready_checked_B_field = np.append(allready_checked_B_field, B_field_name)

    for file in grouped_filenames:
        plot_pair_comparison(file)