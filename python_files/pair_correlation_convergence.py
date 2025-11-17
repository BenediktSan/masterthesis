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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

import header_constants as hc



if os.path.exists("build") == False:
    os.mkdir("build")

if os.path.exists("build/plots") == False:
    os.mkdir("build/plots")

if os.path.exists("build/plots/pair_corr_convergence") == False:
    os.mkdir("build/plots/pair_corr_convergence")

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



def extract_number(s):
    match = re.search(r'\d+', s)
    return int(match.group()) if match else None

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
        if("type"in file_path):
            skip_lines = 2

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



#Functions for 1/N fit

def fit(x, a):
    return a/x 

def calc_fit(x, y):
    # Fit the data to the model function
    param, cov = curve_fit(fit, x, y)

    err = np.sqrt(np.diag(cov))

    return param, err

def calc_difference(final_filename, interim_filename, correlation_type):

    if(extract_number(interim_filename[-7:-1]) != correlation_type):
        print(f"\n\nError: correlation type {correlation_type} does not match file {interim_filename}\n\n")
        return np.array([])

    #get pair correlations

    t_pair = np.array([])

    final_pair_correlations = []

    data_final = read_txt_file("build/daten/"+final_filename)  

    if(len(data_final) == 0):
        return  
    
    #convert to np.array
    final_pair_correlations = []

    for i in range(0, len(data_final)):
        final_pair_correlations.append(data_final[i])

    final_pair_correlations = np.asarray(final_pair_correlations, dtype = np.float32)
    #extract t
    t_pair = final_pair_correlations[-1]
    final_pair_correlations = np.delete(final_pair_correlations,-1,0)



    ##array for correlations
    #set to 0 for all correlations

    interim_pair_correlations = np.array([])

    data_interim = read_txt_file("build/daten/"+interim_filename)  

    if(len(data_interim) == 0):
        return  
    
    #convert to np.array
    interim_pair_correlations =[]
    for i in range(0, len(data_interim)):
        interim_pair_correlations.append(data_interim[i])

    interim_pair_correlations = np.asarray(interim_pair_correlations, dtype = np.float32)
    interim_pair_correlations = np.delete(interim_pair_correlations,-1,0)


    #extract system amount
    system_amount = interim_pair_correlations[0]
    interim_pair_correlations = np.delete(interim_pair_correlations,0,0)
    system_amount = system_amount[system_amount > 0]

    #extract t
    t_pair = interim_pair_correlations[-1]
    pair_correlations = np.delete(interim_pair_correlations,-1,0)



    amount_of_interim_data = len(system_amount)

    if amount_of_interim_data != len(interim_pair_correlations):
        print(f"\n\nError: amount of interim data {amount_of_interim_data} does not match amount of correlations {len(interim_pair_correlations)}\n\n")
        


    deviation = np.zeros(amount_of_interim_data)
    for i in range(0, amount_of_interim_data):
        delta = 0
        step_amount = len(t_pair)
        I = sum(abs(final_pair_correlations[correlation_type]))
        for j in range(0,step_amount):
            delta += abs((interim_pair_correlations[i][j] - final_pair_correlations[correlation_type][j]) / I * 100)
        #maybe even take the square root
        deviation[i] =  delta / step_amount
    return deviation


def plot_correlation_convergence_inlet(pdfname, system_amount, differences, N_fit, fit_params, interim_filenames):

    # Drop-in replacement for your two convergence plotting blocks.
    # Adds a small inset that shows the full data range (overview) and keeps colors aligned.
    # Requires: numpy, matplotlib, mpl_toolkits.axes_grid1.inset_locator.inset_axes



    # --- First plot: differences vs N (with inset showing the whole data range) ---
    size_label = hc.size_label

    fig, ax = plt.subplots()
    plt.rc('axes')

    # plot data and fits, plotting into main axes and inset with matched colors
    ax.axhline(y=0, xmin=0, xmax=1, linestyle="dotted", color="black")  # optional

    axins = inset_axes(ax, width="55%", height="45%", loc='upper right')
    # If you prefer to hide ticks/labels in inset:
    #axins.set_xticks([])
    #axins.set_yticks([])

    # Reset cycles so inset and main use the same color order
    ax.set_prop_cycle(None)
    axins.set_prop_cycle(None)

    # We'll collect linear y-range across all series for the inset limits
    all_y_min = np.inf
    all_y_max = -np.inf

    for i in range(len(differences)):
        # prepare label and x/y arrays
        label = r"$c=$" + str(extract_number(interim_filenames[i][-10:]))

        x = np.asarray(system_amount)
        y = np.asarray(differences[i])

        # plot marker/points on main axes and record the returned Line2D to obtain its color
        main_line, = ax.plot(x, y, "^", label=label)

        # get color from the just-plotted line so the fit uses the same color
        color = main_line.get_color()

        # plot fit curve on main axes using the same color
        fit_line_x = N_fit
        fit_line_y = fit(N_fit, *fit_params[i])
        ax.plot(fit_line_x, fit_line_y, "--", color=color)

       # For the fit: you plotted fit(N_fit, ...) vs 1/N earlier; keep consistent:
        fit_x = 1.0 / np.asarray(N_fit)
        fit_y = fit(N_fit, *fit_params[i])

        # Plot into inset for overview
        x = 1.0/x
        axins.plot(x, y, "^", color=color)
        axins.plot(fit_x, fit_y, "--", color=color)

        # update global y min/max for the inset limits
        all_y_min = min(all_y_min, np.nanmin(y), np.nanmin(fit_y))
        all_y_max = max(all_y_max, np.nanmax(y), np.nanmax(fit_y))

    # tidy up main axes
    ax.set_xlabel(hc.interim_amount_label, fontsize=size_label)
    ax.set_ylabel(hc.deviation_3D_label_pair_corr, fontsize=size_label)
    ax.legend(loc = "center",bbox_to_anchor=(0.85, 0.25), fontsize='small')

    axins.set_xlabel(hc.interim_amount_label_inverse)
    axins.locator_params(axis='x', nbins=7)
    #axins.set_ylabel(r'$\langle \Delta F_x \rangle$ / $\%$')

    # Make the inset show the whole data range (x from main axes and y from aggregated min/max)
    axins.set_xlim(axins.get_xlim())
    # add small padding to avoid data touching the border
    pad = 0.05 * (all_y_max - all_y_min) if (all_y_max > all_y_min) else 0.1
    axins.set_ylim(all_y_min - pad, all_y_max + pad)

    pdfname = pdfname + "_inlet"  # remove "_inverse" suffix for first plot
    outpath = "build/plots/pair_corr_convergence/" + pdfname + ".pdf"
    fig.savefig(outpath, bbox_inches='tight')

    # If you want to replace matched files in the thesis folder
    matching_paths = find_pdfs(pdfname)
    for path in matching_paths:
        print(f"\n\tReplacing plot in master thesis folder: {path}")
        fig.savefig(path, bbox_inches='tight')

    plt.close(fig)

       

def plot_interim_difference_for_pair_corr(filename, interim_results_filename):#,FID_exp, t_exp):


    old_name = []

    
    if "_old" in filename:
        old_name = "_old"
    else:
        old_name = ""



    #get correlations

    ##array forcorrelations
    #set to 0 for all correlations
    t_pair = np.array([])

    correlations = np.array([])

    data = read_txt_file("build/daten/"+interim_results_filename)  

    if(len(data) == 0):
        return  
    
    #convert to np.array
    correlations = np.array([data[0]], dtype=np.float32)
    for i in range(1, len(data)):
        correlations = np.append(correlations,np.array([data[i]],dtype=np.float32), axis = 0)

    #extract system amount
    system_amount = correlations[0]
    correlations = np.delete(correlations,0,0)
    system_amount = system_amount[system_amount > 0]

    #extract t
    t_pair = correlations[-1]
    pair_correlations = np.delete(correlations,-1,0)


    #extract correct pair correlation
    correlation_type = int(extract_number(interim_results_filename[-7:-1]))

    #convert to np.array
    final_pair_correlations = []

    data_final = read_txt_file("build/daten/"+filename)  


    for i in range(0, len(data_final)):
        final_pair_correlations.append(data_final[i])

    final_pair_correlations = np.asarray(final_pair_correlations, dtype = np.float32)
    #extract t
    t_pair = final_pair_correlations[-1]
    final_pair_correlations = np.delete(final_pair_correlations,-1,0)

    #convert t to µs
    t_pair = t_pair/(hbar )*10**-30 * 10**6 

    
    #And now we also do a retrospective correction to the lattice constant, because I used a wrong one in th e simulation
    new_a = 2.72325
    my_a = 2.7315
    factor_a = (new_a/my_a)**(-3)
    #print( f"factor a = {factor_a}")
    t_pair = t_pair / factor_a
    




    size_label = hc.size_label



    plt.figure()
    plt.rc('axes', labelsize=size_label)
    plt.axhline(y = 0, xmin = 0, xmax = 1,linestyle ="dotted", color = "black")
    #converting into µs and adding all factors neglected in simulation
    plt.plot(t_pair , final_pair_correlations[correlation_type], color = hc.FID_color_pair_corr,label=r"Final correlation")
    #plt.plot(t_exp, FID_exp, color ="purple", label=r"experimental data")
    for i in range(0, len(system_amount)):
        plt.plot(t_pair[:len(correlations[i])], correlations[i], label= str(int(system_amount[i])) + " i.c.")
    plt.xlabel( hc.time_label)
    plt.ylabel( hc.pair_correlation_label )
    if(len(system_amount) < 8 ):
        plt.legend()
    pdfname = "pair_corr_comparison_" + interim_results_filename[:-4]
    plt.savefig("build/plots/pair_corr_convergence/" + pdfname + ".pdf")
    #check if plot is in master thesis folder and replace it if it is
    matching_paths = find_pdfs( pdfname)
    for path in matching_paths:
        print(f"\n\tReplacing plot in master thesis folder: {path}")
        plt.savefig(path )
    plt.close()

    return

def plot_pair_correlation_convergence(final_filename, interim_filenames):

    print(f"\nstarted plotting for \n{final_filename}\n\n")




    system_size = extract_number(final_filename[20:])


    pair_correlation_types = []
    correlation_types_name = "corr_types_"


    for i in range(0, len(interim_filenames)):
        pair_correlation_types = np.append(pair_correlation_types, extract_number(interim_filenames[i][-7:-1]))
        correlation_types_name = correlation_types_name + str(int(extract_number(interim_filenames[i][-7:-1]))) + "_"

    pair_correlation_types = np.asarray(pair_correlation_types, dtype = int)



    B_field_names = "[" + final_filename.split("[")[1].split("]")[0] + "]_"

    #calculate differences

    differences = np.array([])
    fit_params = np.array([])
    system_amount = np.array([])


    for i in range(0, len(pair_correlation_types)): 
        file = final_filename
        file_B_field = "[" + file.split("[")[1].split("]")[0] + "]"

        #add B-field
        if(file_B_field not in B_field_names):
            B_field_names = B_field_names + file_B_field + "_"

        if(file[:42] not in interim_filenames[i]):
            print(f"\n\nWarning: File {file} and {interim_filenames[i]} do not match.\n")
            continue
        print(f"\nFile {file} and {interim_filenames[i]} opened\n")

        
        data = read_txt_file("build/daten/"+interim_filenames[i])  

        if(len(data) == 0):
            continue
    
        #convert to np.array
        system_amount = np.array([data[0]], dtype=np.float32)


        system_amount = system_amount[system_amount > 0]
        if(i == 0):
            differences = np.array([calc_difference(file, interim_filenames[i], pair_correlation_types[i])])
        else:
            differences = np.append(differences, np.array([calc_difference(file, interim_filenames[i], pair_correlation_types[i])]), axis = 0)

        #calc fit 
        params, err = calc_fit(system_amount, differences[i])
        print(f"\n\tFit parameters for {str(extract_number(interim_filenames[0][3:]))}sys and {str(extract_number(interim_filenames[0][20:]))}P with f(x) = a/x:")
        print(f"\ta = {params[0]} ± {err[0]}")
        #print(f"\tb = {params[1]} ± {err[1]}")
        if(i== 0):
            fit_params = np.array([params])
        else:
            fit_params = np.append(fit_params, np.array([params]), axis = 0)



    N_fit = np.linspace(system_amount[0], system_amount[-1], 1000)

    size_label = hc.size_label
    plt.figure()
    plt.rc('axes', labelsize=size_label)
    #plt.axhline(y = 0, xmin = 0, xmax = 1,linestyle ="dotted", color = "black")
    for i in range(0, len(differences)):


        file_B_field = "[" + interim_filenames[i].split("[")[1].split("]")[0] + "]"
        if(len(B_field_names) < 10):
            file_B_field  = ""
        plt.plot(system_amount, differences[i ],"^", label=fr"$c = {pair_correlation_types[i]}$")
        #plt.plot(system_amount, differences[i ], color = plt.gca().lines[-1].get_color())
        plt.plot(N_fit, fit(N_fit, *fit_params[i] ),"--", color = plt.gca().lines[-1].get_color())#, label=str(extract_number(filenames_convergence_new[i][7:])) + r"$^3$ Particles fit")
    plt.xlabel( hc.interim_amount_label)
    plt.ylabel( hc.deviation_3D_label_pair_corr )
    plt.legend()    
    pdfname = "interim_pair_deviation_" +str(extract_number(interim_filenames[0][10:])) + "P_" + correlation_types_name + B_field_names + str(extract_number(interim_filenames[0])) +"sys"
    plt.savefig("build/plots/pair_corr_convergence/" + pdfname + ".pdf")
    #check if plot is in master thesis folder and replace it if it is
    matching_paths = find_pdfs( pdfname)
    for path in matching_paths:
        print(f"\n\tReplacing plot in master thesis folder: {path}")
        plt.savefig(path )
    plt.close()

    plot_correlation_convergence_inlet(pdfname, system_amount, differences, N_fit, fit_params, interim_filenames)





all_filenames = []

for file in os.listdir("build/daten"):
    all_filenames.append(file)

all_filenames = np.asarray(all_filenames)

#get rid of pair correlation files



#mask = "pair" not in all_filenames
#all_filenames = all_filenames[mask]
mask = [False] * len(all_filenames)
for i in range(0, len(all_filenames)):
    mask[i] = "pair" in all_filenames[i]
all_filenames = all_filenames[mask]

allready_checked_sys_amount = np.array([])
allready_checked_sys_size = np.array([])
allready_checked_B_field = np.array([])

interim_name = np.array([])
final_name = np.array([])


for file in all_filenames:
    
    system_amount = extract_number(file[:10])
    system_size = extract_number(file[10:])
    #sort for B_fields
    B_field_name = "[" + file.split("[")[1].split("]")[0] + "]"

    continue_flag = True

    for i in range(0, allready_checked_sys_amount.size):
        if system_amount == allready_checked_sys_amount[i] and B_field_name == allready_checked_B_field[i] and system_size == allready_checked_sys_size[i]:
            #print(f"Already checked {system_amount} and {B_field_name}")
            continue_flag = False

    if continue_flag == False:
        continue


    mask = [False] * len(all_filenames)
    for i in range(0, len(all_filenames)):
        if(system_amount == extract_number(all_filenames[i][:10]) and B_field_name in all_filenames[i] and (system_size == extract_number(all_filenames[i][10:])) and "2D" not in all_filenames[i]):
            mask[i] = True
    grouped_filenames = all_filenames[mask]

    #print(f"system amount {system_amount} mask{mask} groups{grouped_filenames}")

    allready_checked_sys_amount = np.append(allready_checked_sys_amount, system_amount)
    allready_checked_B_field = np.append(allready_checked_B_field, B_field_name)
    allready_checked_sys_size = np.append(allready_checked_sys_size, system_size)

    interim_mask = [False] * len(grouped_filenames)
    for i in range(0, len(grouped_filenames)):
        if("_averaged_interim_pair_Correlation_results_type_" in grouped_filenames[i]):
            interim_mask[i] = True

    interim_name = grouped_filenames[interim_mask]

    final_mask = [False] * len(grouped_filenames)
    for i in range(0, len(grouped_filenames)):  
        if("averaged_pair_Correlation" in grouped_filenames[i] and "interim" not in grouped_filenames[i]):
            final_mask[i] = True
            
    final_name = grouped_filenames[final_mask]

    interim_name = np.sort(interim_name)
    final_name = np.sort(final_name)

    if interim_name.size == 0:
        continue

    for correlation_name in final_name:
        system_size = extract_number(correlation_name[20:])
        matching_pair_correlations = interim_name
        mask = [False] * len(interim_name)
        for i in range(0, len(interim_name)):
            if system_size == extract_number(interim_name[i][20:]):
                mask[i] = True
        matching_pair_correlations = matching_pair_correlations[mask]

        print(correlation_name)
        print(matching_pair_correlations)

        plot_pair_correlation_convergence(correlation_name, matching_pair_correlations)
    


    print(final_name)

    for i in range(0, interim_name.size):
        if interim_name.size == 1:
            #plot_interim_difference_for_pair_corr(final_name,interim_name[0])
            continue
        plot_interim_difference_for_pair_corr(final_name[0], interim_name[i])#, exp_111, t_111)                    plot_interim_difference(file, interim_results_file)
