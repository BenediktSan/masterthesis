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
import copy


sys.path.append('/home/sander/Data/Masterarbeit/python_files/')

sys.path.append('/mnt/d/Dateien/Uni/Masterarbeit/python_files/')
import header_constants as hc


if os.path.exists("build") == False:
    os.mkdir("build")

if os.path.exists("build/plots") == False:
    os.mkdir("build/plots")

if os.path.exists("build/plots2D") == False:
    os.mkdir("build/plots2D")

if os.path.exists("build/plots2D/pair_correlation_comparison_external_data") == False:
    os.mkdir("build/plots2D/pair_correlation_comparison_external_data")




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



#Implement taylored correlation functions from stolze89 & Böhm

def get_theoretical_correlation_func(t_pair, corr_type):

    
    allowed_corr_types = np.array([0,1,2,3])    
    if corr_type not in allowed_corr_types:
        if t == 0:
            print("### #### ERROR in calculation of theoretical functions ######")
            print(f"\tcorrelation type {corr_type} not supported")
        return 0
    
    #get only t values before the pair correlation starts to diverge 
    different_cut_offs = np.array([0.9, 0.75, 0.85, 0.85])
    if len(different_cut_offs) != len(allowed_corr_types):
        print("##### ERROR: cut_off lenth array does not match lenth of alowed types")
    t = copy.deepcopy(t_pair)
    mask = t > different_cut_offs[corr_type]
    #t[mask] = 0
    t[mask] = float("NaN")



    #coeff for r = (0,0) correlation
    taylor_coefficents_c0 = np.array([1, 8, 200, 8200, 483160, 39215728, 4278095120, 609534985096])
    #coeff for r = (0,1) correlation
    taylor_coefficents_c1 = np.array([0, -2, -68, -3080, -184000, -14570560 ])
    #coeff for r = (0,2) correlation
    taylor_coefficents_c2 = np.array([0, 0, 6, 410, 28392, 2308476])    
    #coeff for r = (1,1) correlation
    taylor_coefficents_c3 = np.array([0, 0, 12, 760, 48384, 3702888])
 
    all_coefficients = [taylor_coefficents_c0, taylor_coefficents_c1, taylor_coefficents_c2, taylor_coefficents_c3]

    function_value = np.zeros(len(t))

    val = np.zeros(len(t))
    for l in range(0, len(all_coefficients[0])):
        coeff = all_coefficients[0][l]
        val = val +  (-1)**l / np.math.factorial(2 * l) *np.power(t, (2*l))


    
    for l in range(0, len(all_coefficients[corr_type])):
        coeff = all_coefficients[corr_type][l]
        function_value = function_value +  (-1)**l / np.math.factorial(2 * l) * coeff * np.power(t, (2*l))

    return function_value


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

def find_data_for_comparison_przemek(filename):

    B_field_name = "[" + filename.split("[")[1].split("]")[0] + "]"

    #find Przemyslaws data

    if(B_field_name == "[001]"):
        new_B_field_name = "100"
    elif(B_field_name == "[011]"):
        print(f"\n##############ERROR FOR THIS B_FIELD {B_field_name} NO DATA EXISTS")
    elif(B_field_name == "[111]"):
        print(f"\n##############ERROR FOR THIS B_FIELD {B_field_name} NO DATA EXISTS")

    przemek_t = []
    przemek_auto_correlation = []
    przemek_pair_correlations = []

    #add buffer because site =6 
    pair_corr_buffer = []
    spin_site = []


    for file in os.listdir("../theoretical_data/2D_Heisenberg/Przemek/"):
        if ".hdf5" in file:
            przemek_data = h5py.File("../theoretical_data/2D_Heisenberg/Przemek/"+file, "r")

            przemek_corr = przemek_data['results']['Re_correlation'][0]
            T_max = przemek_data['parameters'].attrs['Tmax']
            time_points = przemek_data['parameters'].attrs['num_TimePoints']

            t = np.linspace(0, T_max, time_points)
            
            spin_site_num = przemek_data['parameters'].attrs['spin_site']

            przemek_t.append(t)

            #print(spin_site_num)
            
            if spin_site_num == 0:
                przemek_auto_correlation = przemek_corr
            else:
                przemek_pair_correlations.append(przemek_corr)


    przemek_pair_correlations = np.asarray(przemek_pair_correlations)
    przemek_auto_correlation = np.asarray(przemek_auto_correlation)
    przemek_t = np.asarray(przemek_t)

    return przemek_t, przemek_auto_correlation, przemek_pair_correlations

def find_data_for_comparison_timo(filename):

    B_field_name = "[" + filename.split("[")[1].split("]")[0] + "]"

    #find Przemyslaws data

    if(B_field_name == "[001]"):
        new_B_field_name = "100"
    elif(B_field_name == "[011]"):
        print(f"\n##############ERROR FOR THIS B_FIELD {B_field_name} NO DATA EXISTS")
    elif(B_field_name == "[111]"):
        print(f"\n##############ERROR FOR THIS B_FIELD {B_field_name} NO DATA EXISTS")

    timo_t = []
    timo_auto_correlation = []
    timo_pair_correlations = []

    #add buffer because site =6 
    pair_corr_buffer = []
    spin_site = []


    for file in os.listdir("../theoretical_data/2D_Heisenberg/Timo/data/"):
        if ".hdf5" in file:

            timo_data = h5py.File("../theoretical_data/2D_Heisenberg/Timo/data/"+file, "r")

            spin_site_num = extract_number(timo_data['parameters'].attrs['config_file'])

            if spin_site_num == 0:
                timo_corr = timo_data['results']['correlation']["1-1"][0]
                timo_auto_correlation = timo_corr
            else:
                timo_corr = timo_data['results']['correlation']["1-2"][0]
                timo_pair_correlations.append(timo_corr)



            delta_t = timo_data['parameters'].attrs['delta_t']
            time_points = timo_data['parameters'].attrs['num_TimePoints']

            t = np.linspace(0, delta_t*time_points, time_points)


            timo_t.append(t)


    timo_pair_correlations = np.asarray(timo_pair_correlations)
    timo_auto_correlation = np.asarray(timo_auto_correlation)
    timo_t = np.asarray(timo_t)

    return timo_t, timo_auto_correlation, timo_pair_correlations



def plot_pair_comparison(pair_filename):
 

    print(f"\nstarted plotting for \n{pair_filename}\n\n")


    B_field_name =  "[" + pair_filename.split("[")[1].split("]")[0] + "]"


    system_size = extract_number(pair_filename[18:38])
    system_size_name = str(system_size) + r"P"
    print(f"\t\tSystem Size: {system_size} Particles")

    #get pair correlations

    ##array for pair correlations
    #set to 0 for all correlations
    #this cutoff defines the amunt of correlations per subplot
    #this excludes the auto correlation plot
    correlation_cut_off = 4
    sum_of_all_correlations =  np.array([])
    t_pair = np.array([])

    pair_correlations = np.array([])

    data = read_txt_file("build/daten_2D/"+pair_filename)  

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

    auto_correlation = pair_correlations[0]
    pair_correlations = np.delete(pair_correlations,0,0)


    if (len(pair_correlations) < 4):
        print(f"\n\t\tWARNING: Not enough pair correlations found for {pair_filename} to plot all subplots")
        print(f"\n\t\tAdding {(4) - len(pair_correlations)} empty pair correlation  to fill up subplots")
        for i in range(len(pair_correlations), 4):
            pair_correlations = np.append(pair_correlations, np.zeros((1, len(t_pair))), axis=0)

    ##normalizing
    
    norm =  sum_of_all_correlations[0]

    #norm = 5 ** 3 * 8
    sum_of_all_correlations = sum_of_all_correlations / norm
    auto_correlation = auto_correlation / norm
    pair_correlations = pair_correlations / norm

    ##get rid of multiplicity in my data
    pair_correlations /= 4



    ##extract przemeks data from the hdf5 file

    przemek_t, przemek_auto, przemek_pair = find_data_for_comparison_przemek(pair_filename)

    ##normalize przemeks data
    norm =przemek_auto[0]
    przemek_auto = przemek_auto / norm
    przemek_pair = przemek_pair / norm
    #for i in range(0, len(przemek_pair)):
        #przemek_pair[i] = przemek_pair[i] 


    ##extract timos data
    timo_t, timo_auto, timo_pair = find_data_for_comparison_timo(pair_filename)

    ##normalize timos data
    norm =timo_auto[0]
    timo_auto = timo_auto / norm
    timo_pair = timo_pair / norm

    ##rescaling time because of false hamiltonian
    przemek_t /=4
    t_theo = copy.deepcopy(t_pair)
    t_pair /=2
    timo_t /=4

    ##check which pair correlations need to be plotted

    #????????????????
    pair_correlation_types = np.array([0,1,2,3])



    #print(f"\t\tprzemek data has {len(przemek_pair)} pair correlations and {len(przemek_t)} time steps")
    size_label = hc.size_label

    #print(f"\t\tsize of t {t.size} and size of pair_corrs {pair_correlations[0].size} {pair_correlations[1].size}")
    pastel1 = plt.get_cmap('Pastel1')
    pastel2 = plt.get_cmap('Pastel2')
    pastel1.colors += pastel2.colors


    #####NOTE: the [1] correlation calculated by me is (1,1) and [2] is (0,2) because of the way i order my correlations
    #because of this they need to be switched to match the order of the theoretical data
    #I assume the same to be true for timos data




    fig, axs = plt.subplots(2, 2, sharex=True, figsize=(10, 8))
    #fig.suptitle(r"Pair Correlation Comparison " + str(extract_number(pair_filename[10:])) + r" i.c. and " +str(extract_number(pair_filename[20:])) +  r"$^2$ lattice sites")
    axs[0, 0].plot(t_pair, auto_correlation, color = hc.pair_color_2D, label =r"Semi-classical $G_0^{xx}$")
    axs[0, 0].plot(przemek_t[0], przemek_auto, color = hc.przemek_color, label = r"Chebyshec expansion")
    axs[0, 0].plot(timo_t[0], timo_auto, color = hc.timo_color_FID, label = r"nl-spinDMFT")

    axs[0, 0].plot(t_theo, get_theoretical_correlation_func(t_theo, corr_type=0), linestyle = "--", color = hc.time_exp_color, label = r"Time expansion")

    axs[0,0].legend()
    axs[0, 0].set_title('Normalized auto correlation')

    c = pair_correlation_types[1]
    axs[0, 1].plot(t_pair, pair_correlations[0 ], color = hc.pair_color_2D, label = r"Semi-classical $G_" + f"{c}" + r"^{xx}$")
    axs[0, 1].plot(przemek_t[1], przemek_pair[0 ], color = hc.przemek_color, label = r"Chebyshec expansion")
    axs[0, 1].plot(timo_t[1],timo_pair[0 ], color = hc.timo_color_FID, label = r"nl-spinDMFT")
    axs[0, 1].plot(t_theo, get_theoretical_correlation_func(t_theo, corr_type=pair_correlation_types[1]), linestyle = "--", color = hc.time_exp_color, label = r"Time expansion " )
    #axs[0,1].set_xscale("log")
    #axs[0,1].set_yscale("log")
    axs[0, 1].legend()
    axs[0,1].set_title(r"Normalized pair correlation for $\vec{v} = (0,1)$")

        
    c = pair_correlation_types[2]
    axs[1, 0].plot(t_pair,pair_correlations[2], color = hc.FID_color, label = r"Semi-classical $ G_{" + f"{c}" + r"}^{xx}$")
    axs[1, 0].plot(przemek_t[2],przemek_pair[1], color = hc.przemek_color, label = r"Chebyshec expansion")
    axs[1, 0].plot(timo_t[2], timo_pair[2 ], color = hc.timo_color_FID, label = r"nl-spinDMFT")
    axs[1, 0].plot(t_theo, get_theoretical_correlation_func(t_theo, corr_type=pair_correlation_types[2]), color = hc.time_exp_color, linestyle = "--", label = r"Time expansion " )
    #axs[1,0].set_xscale("log")
    #axs[1,0].set_yscale("log")
    axs[1,0].legend()
    axs[1,0].set_title(r'Normalized pair correlation for $\vec{v} = (0,2)$')

    
    c = pair_correlation_types[3]
    axs[1, 1].plot(t_pair, pair_correlations[1], color = hc.pair_color_2D, label = r"Semi-classical $ G_{" + f"{c}" + r"}^{xx}$")
    #axs[1, 1].plot(t_pair, pair_correlations[0], label = r"classical $ G_{" + f"{c}" + r"}^{xx}$")
    #axs[1, 1].plot(t_pair, pair_correlations[2], label = r"classical $ G_{" + f"{c}" + r"}^{xx}$")
    axs[1, 1].plot(przemek_t[3], przemek_pair[2], color = hc.przemek_color, label = r"Chebyshec expansion")
    axs[1, 1].plot(timo_t[3], timo_pair[1 ], color = hc.timo_color_FID, label = r"nl-spinDMFT")
    axs[1, 1].plot(t_theo, get_theoretical_correlation_func(t_theo, corr_type=3), linestyle = "--", color = hc.time_exp_color, label = r"Time expansion" )
    #axs[1,1].set_xscale("log")
    #axs[1,1].set_yscale("log")
    axs[1, 1].legend()
    axs[1,1].set_title(r'Normalized pair correlation for $\vec{v} = (1,1)$')

    #for ax in axs.flat:
    #    ax.set(xlabel=r'$t$ / $J$', ylabel=r'$FID$ ')

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    #for ax in axs.flat:
    #    ax.label_outer()

    axs[1, 0].set_xlabel(hc.time_label_2D,fontsize=size_label)
    axs[1, 1].set_xlabel(hc.time_label_2D,fontsize=size_label)
    axs[0, 0].set_ylabel(hc.certain_pair_correlation_label,fontsize=size_label)
    axs[1, 0].set_ylabel(hc.certain_pair_correlation_label,fontsize=size_label)

    fig.tight_layout()

    fig.savefig("build/plots2D/pair_correlation_comparison_external_data/pair_correlation_2D_comparison_external_data_" + system_size_name + "_" + str(extract_number(pair_filename[10:])) +"sys_" + B_field_name + ".pdf")
    pdfname = "pair_correlation_2D_comparison_external_data_" + system_size_name + "_" + str(extract_number(pair_filename[10:])) +"sys_" + B_field_name 
    plt.savefig("build/plots2D/pair_correlation_comparison_external_data/" + pdfname + ".pdf")
    #check if plot is in master thesis folder and replace it if it is
    matching_paths = find_pdfs( pdfname)
    for path in matching_paths:
        print(f"\n\tReplacing plot in master thesis folder: {path}")
        plt.savefig(path )    
    plt.close()
    return




all_filenames = np.array([])

for file in os.listdir("build/daten_2D"):
    all_filenames = np.append(all_filenames,file)



mask = [False] * len(all_filenames)
for i in range(0, len(all_filenames)):
    mask[i] = "pair" in all_filenames[i] and "2D" in all_filenames[i] and "interim" not in all_filenames[i] and "variance" not in all_filenames[i] and "XXZ" not in all_filenames[i] and "dipolar" not in all_filenames[i]
all_filenames = all_filenames[mask]

# just so its not empty and i can append
allready_checked_sys_amount = np.array([])
allready_checked_B_field = np.array([])


for file in all_filenames:
    
    system_amount = extract_number(file[10:])

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
        if(system_amount == extract_number(all_filenames[i][10:]) and B_field_name in all_filenames[i]):
            mask[i] = True
    grouped_filenames = all_filenames[mask]

    #print(f"system amount {system_amount} mask{mask} groups{grouped_filenames}")

    allready_checked_sys_amount = np.append(allready_checked_sys_amount, system_amount)
    allready_checked_B_field = np.append(allready_checked_B_field, B_field_name)

    print(grouped_filenames)
    for file in grouped_filenames:
        plot_pair_comparison(file)