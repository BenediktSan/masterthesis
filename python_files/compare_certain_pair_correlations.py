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



if os.path.exists("build/plots/pair_comparison") == False:
    os.mkdir("build/plots/pair_comparison")


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


def find_data_for_comparison(filename):

    B_field_name = "[" + filename.split("[")[1].split("]")[0] + "]"

    #find Gräßers data

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


def plot_certain_pair_comparison(grouped_filenames):


    dimension_name = []

    grouped_filenames = sorted(grouped_filenames, key=lambda x: extract_number(x[18:38]) )
    grouped_filenames = np.asarray(grouped_filenames)
 
    print(f"\nstarted plotting for \n{grouped_filenames}\n\n")


    B_field_name =  "[" + grouped_filenames[0].split("[")[1].split("]")[0] + "]"


    system_size = []
    system_size_name = []
    all_system_sizes_names = ""

    for i in range(0, len(grouped_filenames)):
        system_size = np.append(system_size, extract_number(grouped_filenames[i][18:38]))
        system_size_name = np.append(system_size_name, str(int(system_size[i])) )
        all_system_sizes_names += system_size_name[i] + "_"

    sytsem_size = np.asarray(system_size)
    system_size_name = np.asarray(system_size_name)



    
    #get pair correlations

    ##array for pair correlations
    #set to 0 for all correlations
    #this cutoff defines the amunt of correlations per subplot
    #this excludes the auto correlation plot
    t_pair = []

    all_pair_correlations = []

    for i in range(0, len(grouped_filenames)):
    
        data = read_txt_file("build/daten/"+grouped_filenames[i]  )
        if(len(data) == 0):
            return
    

        
        pair_correlations_temp = []
        for j in range(0, len(data)):
            pair_correlations_temp.append(np.array(data[j], dtype=np.float32))

        #ectract time values
        t_pair.append(pair_correlations_temp[-1])
        pair_correlations_temp = np.delete(pair_correlations_temp, -1, 0)
    
        pair_correlations_temp = np.asarray(pair_correlations_temp)
    
        all_pair_correlations.append(pair_correlations_temp)



    t_pair = np.asarray(t_pair, dtype=object)

    #convert t to µs
    for i in range(0, len(t_pair)):
        t_pair[i] = t_pair[i]/(hbar )*10**-30 * 10**6 
    
    #And now we also do a retrospective correction to the lattice constant, because I used a wrong one in th e simulation
    new_a = 2.72325
    my_a = 2.7315
    factor_a = (new_a/my_a)**(-3)
    #print( f"factor a = {factor_a}")
    for i in range(0, len(t_pair)):
        t_pair[i] = t_pair[i] / factor_a
    

    chosen_correlations = np.array([0,1,2,3])

    if len(chosen_correlations) != 4:
        print(f"\n\t ERROR: Code is not written to account for more or less than 4 subplots. Chosen correlations are {chosen_correlations}\n")

    for i in range(0, grouped_filenames.size):
        if (len(all_pair_correlations[i]) < len(chosen_correlations)):

            print(f"\n\t\tWARNING: Not enough pair correlations found for {grouped_filenames[i]} to plot all subplots")
            print(f"\n\t\tAdding {(len(chosen_correlations)) - len(all_pair_correlations[i])} empty pair correlation  to fill up subplots")

            for j in range(len(all_pair_correlations[i]),len(chosen_correlations)):
                all_pair_correlations[i] = np.append(all_pair_correlations[i], np.zeros((1, len(t_pair[i]))), axis=0)

    ##normalizing
    
    normalization_constants = np.zeros(len(grouped_filenames))

    all_pair_correlations = np.asarray(all_pair_correlations, dtype = object)


    for i in range(0, len(grouped_filenames)):
        for j in range(0, len(all_pair_correlations[i])):
            normalization_constants[i] += all_pair_correlations[i][j][0]

    print(f"Normalization constants are {normalization_constants}")

    for i in range(0, len(grouped_filenames)):
        all_pair_correlations[i] /= normalization_constants[i]


    #norm = 5 ** 3 * 8
    #sum_of_all_correlations = sum_of_all_correlations / norm
    #auto_correlation = auto_correlation / norm
    #pair_correlations = pair_correlations / norm



    ##get timos data to also compare it to that
    timo_t, timo_auto, timo_pair = find_data_for_comparison(grouped_filenames[0])


    #ugly workaround for ugly sorting in c++ code
    #for i in range(0, len(system_size)):
    #    if (system_size[i] > 6.):
    #        buffer = all_pair_correlations[i][3]
    #        all_pair_correlations[i][3] = all_pair_correlations[i][4]
    #        all_pair_correlations[i][4] = buffer


    size_label = hc.size_label

    
    fig, axs = plt.subplots(2, 2, sharex=True, figsize=(10, 8))


    #fig.suptitle(r"Pair Correlation Comparison " + str(extract_number(pair_filename[10:])) + r" i.c. and " +str(extract_number(pair_filename[20:])) +  r"$^2$ Particles")
    for i in range(0, len(grouped_filenames)):
        print(grouped_filenames[i])
        axs[0, 0].plot(t_pair[i], all_pair_correlations[i][int(chosen_correlations[0])], label = f"{system_size_name[i]}" + r"$^{3}$ lattice sites")
    axs[0, 0].plot(timo_t,timo_auto, color = hc.timo_color_FID, linestyle = "--", label = r"nl-spinDMFT $G_0^{xx}$")
    #axs[0, 0].plot(timo_t,timo_auto, color = axs[0, 0].lines[-1].get_color(), linestyle = "--", label = r"$m_0 G_0^{xx}$")
    axs[0,0].legend()
    axs[0, 0].set_title('Auto correlation')
    for i in range(0, len(grouped_filenames)):
        axs[0, 1].plot(t_pair[i], all_pair_correlations[i][chosen_correlations[1]], label =  f"{system_size_name[i]}" + r"$^{3}$ lattice sites")
    axs[0, 1].plot(timo_t,timo_pair[chosen_correlations[1] - 1], color = hc.timo_color_FID, linestyle = "--", label = r"nl-spinDMFT " + r"$m_{" + f"{chosen_correlations[1]}"  + r"} G_{" + f"{chosen_correlations[1]}" + r"}^{xx}$")

    #axs[0, 0].plot(timo_t,timo_auto, color = axs[0, 0].lines[-1].get_color(), linestyle = "--", label = r"$m_0 G_0^{xx}$")
    axs[0,1].legend()
    axs[0,1].set_title(r'Pair correlation for $c = $' + str(chosen_correlations[1]))

    for i in range(0, len(grouped_filenames)):
        axs[1, 0].plot(t_pair[i], all_pair_correlations[i][chosen_correlations[2]], label =  f"{system_size_name[i]}" + r"$^{3}$ lattice sites")
    axs[1, 0].plot(timo_t,timo_pair[chosen_correlations[2] - 1], color = hc.timo_color_FID, linestyle = "--", label = r"nl-spinDMFT " + r"$m_{" + f"{chosen_correlations[2]}"  + r"} G_{" + f"{chosen_correlations[2]}" + r"}^{xx}$")
    #axs[0, 0].plot(timo_t,timo_auto, color = axs[0, 0].lines[-1].get_color(), linestyle = "--", label = r"$m_0 G_0^{xx}$")
    axs[1,0].legend()
    axs[1,0].set_title(r'Pair correlation for $c = $' + str(chosen_correlations[2]))
    for i in range(0, len(grouped_filenames)):
        axs[1, 1].plot(t_pair[i], all_pair_correlations[i][chosen_correlations[3]], label =  f"{system_size_name[i]}" + r"$^{3}$ lattice sites")
        #axs[1, 1].plot(t_pair[i], all_pair_correlations[i][4], label = r"Pair correlation for " + f"{system_size_name[i]}" + r"$^{3}$")
        #axs[1, 1].plot(t_pair[2], all_pair_correlations[2][chosen_correlations[3]], label = r"Pair correlation for " + f"{system_size_name[2]}" + r"^{2}$")
    axs[1, 1].plot(timo_t,timo_pair[chosen_correlations[3] - 1], color = hc.timo_color_FID, linestyle = "--", label = r"nl-spinDMFT " + r"$m_{" + f"{chosen_correlations[3]}"  + r"} G_{" + f"{chosen_correlations[3]}" + r"}^{xx}$")
    #axs[1, 1].plot(timo_t,timo_pair[3], color = axs[1, 0].lines[-1].get_color(), linestyle = "--", label = r"Gräßers $G_0^{xx}$")

    axs[1,1].legend()
    axs[1,1].set_title(r'Pair correlation for $c = $' + str(chosen_correlations[3]))

    #axs[1, 1].set_title('C')

    #for ax in axs.flat:
    #    ax.set(xlabel=r'$t$ / $J$', ylabel=r'$FID$ ')
    axs[1, 0].set_xlabel(hc.time_label,fontsize=size_label)
    axs[1, 1].set_xlabel(hc.time_label,fontsize=size_label)
    axs[0, 0].set_ylabel(hc.pair_correlation_label,fontsize=size_label)
    axs[1, 0].set_ylabel(hc.pair_correlation_label,fontsize=size_label)

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    #for ax in axs.flat:
    #    ax.label_outer()

    fig.tight_layout()
    pdfname = all_system_sizes_names + str(extract_number(grouped_filenames[0][:10])) +"sys_" + B_field_name
    plt.savefig("build/plots/pair_comparison/compare_certain_pair_corrs_" + pdfname + ".pdf")
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
    mask[i] = "pair" in all_filenames[i] and "2D" not in all_filenames[i] and "interim" not in all_filenames[i] and "variance" not in all_filenames[i]
all_filenames = all_filenames[mask]

# just so its not empty and i can append
allready_checked_sys_amount = np.array([])
allready_checked_B_field = np.array([])


for file in all_filenames:
    
    system_amount = extract_number(file[:10])

    continue_flag = True

    #sort for B_fields
    B_field_name = "[" + file.split("[")[1].split("]")[0] + "]"

    for i in range(0, allready_checked_sys_amount.size):
        if system_amount == allready_checked_sys_amount[i] and B_field_name == allready_checked_B_field[i]:
            #print(f"Already checked {system_amount} and {B_field_name}")
            continue_flag = False

    if continue_flag == False:
        continue

    
    mask = [False] * len(all_filenames)
    for i in range(0, len(all_filenames)):
        if(system_amount == extract_number(all_filenames[i][:10]) and B_field_name in all_filenames[i]):
            mask[i] = True
    grouped_filenames = all_filenames[mask]


    #print(f"system amount {system_amount} mask{mask} groups{grouped_filenames}")

    allready_checked_sys_amount = np.append(allready_checked_sys_amount, system_amount)
    allready_checked_B_field = np.append(allready_checked_B_field, B_field_name)

    if(len(grouped_filenames ) < 1):
        continue

    plot_certain_pair_comparison(grouped_filenames)