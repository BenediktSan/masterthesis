from math import sqrt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
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

if os.path.exists("build/plots/convergence") == False:
    os.mkdir("build/plots/convergence")

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

def fit_sqrt(x, a):
    return a/np.sqrt(x)

def calc_fit(x, y):
    # Fit the data to the model function
    param, cov = curve_fit(fit, x, y)

    err = np.sqrt(np.diag(cov))

    return param, err

def calc_fit_sqrt(x, y):
    # Fit the data to the model function
    param, cov = curve_fit(fit_sqrt, x, y)

    err = np.sqrt(np.diag(cov))

    return param, err

def calc_difference(final_filename, interim_filename):

    #get correlations

    final_correlations, err, t= np.genfromtxt("build/daten/"+final_filename,skip_header= 1, unpack = True)

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
    

    if("_old" in final_filename):
        t_tilde /= 2    

    ##array forcorrelations
    #set to 0 for all correlations
    t_pair = np.array([])

    correlations = np.array([])

    data = read_txt_file("build/daten/"+interim_filename)  

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



    amount_of_files = len(system_amount)



    deviation = np.zeros(amount_of_files)
    for i in range(0, amount_of_files):
        delta = 0
        step_amount = len(t_tilde)
        I = sum(abs(final_correlations))
        for j in range(0,step_amount):
            delta += abs((correlations[i][j] - final_correlations[j]) / I * 100)
        #maybe even take the square root
        deviation[i] =  delta / step_amount
    return deviation


def plot_correlation_convergence_inlet(pdfname, system_amount, differences, N_fit, fit_params, filenames_convergence_new, B_field_names):

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
        file_B_field = "[" + filenames_convergence_new[i].split("[")[1].split("]")[0] + "]"
        if len(B_field_names) < 10:
            file_B_field = ""

        B_field_label = ""
        if file_B_field == "[001]":
            B_field_label = hc.B001
        elif file_B_field == "[111]":
            B_field_label = hc.B111
        elif file_B_field == "[011]":
            B_field_label = hc.B011
        else:
            print(f"Warning: Unknown B-field direction {B_field_names}.\n")



        #label = (str(extract_number(filenames_convergence_new[i][7:])) + r"$^3$ lattice sites with " + B_field_label)
        label =  B_field_label

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
    ax.set_xlabel( hc.interim_amount_label, fontsize=hc.size_label)
    ax.set_ylabel( hc.deviation_3D_label_interim, fontsize=hc.size_label)
    ax.legend(loc = "center",bbox_to_anchor=(0.85, 0.25), fontsize='small')

    axins.set_xlabel(hc.interim_amount_label_inverse)
    axins.locator_params(axis='x', nbins=7)
    #axins.set_ylabel(r'$\langle \Delta F_x \rangle$ / $\%$')

    # Make the inset show the whole data range (x from main axes and y from aggregated min/max)
    axins.set_xlim(axins.get_xlim())
    # add small padding to avoid data touching the border
    pad = 0.05 * (all_y_max - all_y_min) if (all_y_max > all_y_min) else 0.1
    axins.set_ylim(all_y_min - pad, all_y_max + pad)

    pdfname = pdfname[:-8] + "_inlet"  # remove "_inverse" suffix for first plot
    outpath = "build/plots/convergence/" + pdfname + ".pdf"
    fig.savefig(outpath, bbox_inches='tight')

    # If you want to replace matched files in the thesis folder
    matching_paths = find_pdfs(pdfname)
    for path in matching_paths:
        print(f"\n\tReplacing plot in master thesis folder: {path}")
        fig.savefig(path, bbox_inches='tight')

    plt.close(fig)

       


def plot_interim_difference(filename, interim_results_filename):#,FID_exp, t_exp):




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


    
    C, err, t= np.genfromtxt("build/daten/"+filename,skip_header= 1, unpack = True)

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
    #print(f"\n\t\tstepsize for t_tilde = {t_tilde[1]}")

    if("_old" in filename):
        t_tilde /= 2


    size_label = hc.size_label

    #print(f"\t\tsize of t_tilde {t_tilde.size} and size of pair_corrs {pair_correlations[0].size} {pair_correlations[1].size}")


    plt.figure()
    plt.rc('axes', labelsize=size_label)
    plt.axhline(y = 0, xmin = 0, xmax = 1,linestyle ="dotted", color = "black")
    #converting into µs and adding all factors neglected in simulation
    plt.plot(t_tilde , C, color = "green",label=r"Correlation")
    plt.fill_between(t_tilde, C + err, C - err, color = "green", alpha = 0.3, label ="estimated averaging error")
    #plt.plot(t_tilde_2 , C, label=r"C3", alpha = 0.5)
    #plt.plot(t_exp, FID_exp, color ="purple", label=r"experimental data")
    for i in range(0, len(system_amount)):
        plt.plot(t_tilde[:len(correlations[i])], correlations[i], label= str(int(system_amount[i])) + " i.c.")
    plt.xlabel( hc.time_label)
    plt.ylabel( hc.FID_label )
    if(len(system_amount) < 8 ):
        plt.legend()
    plt.tight_layout()
    pdfname = "comparison_" + interim_results_filename[:-4] 
    plt.savefig("build/plots/convergence/" + pdfname + ".pdf")
    #check if plot is in master thesis folder and replace it if it is
    matching_paths = find_pdfs( pdfname)
    for path in matching_paths:
        print(f"\n\tReplacing plot in master thesis folder: {path}")
        plt.savefig(path )
    plt.close()

    return

def plot_correlation_convergence(final_filenames, interim_filenames):

    print(f"\nstarted plotting for \n{final_filenames}\n\n")
    if len(final_filenames) != len(interim_filenames):
        print(f"Warning: Different amount of final and interim files.\n")
        print(f"Final files: {final_filenames}\nInterim files: {interim_filenames}\n")


    old_name = []

    for file in final_filenames:
        if "_old" in file:
            old_name.append("old")
        else:
            old_name.append("")

    old_name = np.asarray(old_name)



    
    #create filenames

    filenames_convergence_new = np.array([])
    for i in range(0, len(interim_filenames)):
        filenames_convergence_new = np.append(filenames_convergence_new, interim_filenames[i])


    if("lido3" in str(filenames_convergence_new)):
        for i in range(0,filenames_convergence_new.size):
            filenames_convergence_new[i] = filenames_convergence_new[i].replace("lido3_", "")


    B_field_names = "[" + final_filenames[0].split("[")[1].split("]")[0] + "]_"
    B_field_label = ""
    if B_field_names == "[001]_":
        B_field_label = hc.B001
    elif B_field_names == "[111]_":
        B_field_label = hc.B111
    elif B_field_names == "[011]_":
        B_field_label = hc.B011
    else:
        print(f"Warning: Unknown B-field direction {B_field_names}.\n")
    #calculate differences

    differences =[]
    fit_params = np.array([])
    system_amount = np.array([])


    for i in range(0, len(final_filenames)): 
        file = final_filenames[i]
        file_B_field = "[" + file.split("[")[1].split("]")[0] + "]"
        system_size_temp = extract_number(file[10:])



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
        differences.append(calc_difference(file, interim_filenames[i]))


        #calc fit 
        
        params, err = calc_fit(system_amount, differences[i])
        print(f"\n\tFit parameters for {str(extract_number(interim_filenames[0][:-7]))}sys and {str(extract_number(filenames_convergence_new[0][8:]))}P with f(x) = a/x:")
        print(f"\ta = {params[0]} ± {err[0]}")
        #print(f"\tb = {params[1]} ± {err[1]}")
        if(i== 0):
            fit_params = np.array([params])
        else:
            fit_params = np.append(fit_params, np.array([params]), axis = 0)

    differences = np.asarray(differences)

    N_fit = np.linspace(system_amount[0], system_amount[-1], 1000)

    size_label = hc.size_label
    plt.figure()
    plt.rc('axes', labelsize=size_label)
    #plt.axhline(y = 0, xmin = 0, xmax = 1,linestyle ="dotted", color = "black")
    for i in range(0, len(differences)):

        file_B_field = "[" + filenames_convergence_new[i].split("[")[1].split("]")[0] + "]"
        if(len(B_field_names) < 10):
            file_B_field  = ""
        plt.plot(system_amount, differences[i ],"^", label=str(extract_number(filenames_convergence_new[i][7:])) + r"$^3$ lattice sites for " + B_field_label + old_name[i])
        #plt.plot(system_amount, differences[i ], color = plt.gca().lines[-1].get_color())
        plt.plot(N_fit, fit(N_fit, *fit_params[i] ),"--", color = plt.gca().lines[-1].get_color())#, label=str(extract_number(filenames_convergence_new[i][7:])) + r"$^3$ Particles fit")
    plt.xlabel( hc.interim_amount_label)
    plt.ylabel( hc.deviation_3D_label_interim )
    plt.legend()    
    plt.tight_layout()

    pdfname = "interim_deviation_" + B_field_names + str(extract_number(interim_filenames[0][:-7]))
    plt.savefig("build/plots/convergence/" + pdfname + ".pdf")
    #check if plot is in master thesis folder and replace it if it is
    matching_paths = find_pdfs( pdfname)
    for path in matching_paths:
        print(f"\n\tReplacing plot in master thesis folder: {path}")
        plt.savefig(path )
    plt.close()





    plt.figure()
    plt.rc('axes', labelsize=size_label)
    #plt.axhline(y = 0, xmin = 0, xmax = 1,linestyle ="dotted", color = "black")
    for i in range(0, len(differences)):

        file_B_field = "[" + filenames_convergence_new[i].split("[")[1].split("]")[0] + "]"
        if(len(B_field_names) < 10):
            file_B_field  = ""
        plt.plot(1/(system_amount ), differences[i ],"^", label=str(extract_number(filenames_convergence_new[i][7:])) + r"$^3$ Particles " + B_field_label + old_name[i])
        #plt.plot(system_amount, differences[i ], color = plt.gca().lines[-1].get_color())
        plt.plot(1/(N_fit), fit(N_fit, *fit_params[i] ),"--", color = plt.gca().lines[-1].get_color())#, label=str(extract_number(filenames_convergence_new[i][7:])) + r"$^3$ Particles fit")
    plt.xlabel( hc.interim_amount_label_inverse)
    plt.locator_params(axis='x', nbins=7)    #plt.xscale('function', functions=(inverse, inverse))
    plt.ylabel( hc.deviation_3D_label_interim)
    plt.legend()    
    plt.tight_layout()

    pdfname = "interim_deviation_" + B_field_names + str(extract_number(interim_filenames[0][:-7])) + "_inverse"
    plt.savefig("build/plots/convergence/" + pdfname + ".pdf")
    #check if plot is in master thesis folder and replace it if it is
    matching_paths = find_pdfs( pdfname)
    for path in matching_paths:
        print(f"\n\tReplacing plot in master thesis folder: {path}")
        plt.savefig(path )
    plt.close()


    plot_correlation_convergence_inlet(pdfname, system_amount, differences, N_fit, fit_params, filenames_convergence_new, B_field_names)

    return

def plot_correlation_convergence_sqrt(final_filenames,interim_filenames):

    if("no_time_int" not in str(final_filenames)):
        return

    print(f"\nstarted plotting for \n{final_filenames}\n\n")



    
    #create filenames

    filenames_convergence_new = np.array([])
    for i in range(0, len(interim_filenames)):
        filenames_convergence_new = np.append(filenames_convergence_new, interim_filenames[i])




    B_field_names = "[" + final_filenames[0].split("[")[1].split("]")[0] + "]_"
    B_field_label = ""
    if B_field_names == "[001]_":
        B_field_label = hc.B001
    elif B_field_names == "[111]_":
        B_field_label = hc.B111
    elif B_field_names == "[011]_":
        B_field_label = hc.B011
    else:
        print(f"Warning: Unknown B-field direction {B_field_names}.\n")
    #calculate differences
    pdfname = "interim_deviation_" + B_field_names + str(extract_number(interim_filenames[0][:-7])) + "_no_time_average"

    differences =[]
    fit_params = np.array([])
    system_amount = np.array([])


    for i in range(0, len(final_filenames)): 
        file = final_filenames[i]
        file_B_field = "[" + file.split("[")[1].split("]")[0] + "]"
        system_size_temp = extract_number(file[10:])



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
        differences.append(calc_difference(file, interim_filenames[i]))


        #calc fit 
        
    differences = np.asarray(differences)

    N_fit = np.linspace(system_amount[0], system_amount[-1], 1000)


    params, err = calc_fit_sqrt(system_amount, differences[i])
    print(f"\n\tFit parameters for {str(extract_number(interim_filenames[0][:-7]))}sys and {str(extract_number(filenames_convergence_new[0][8:]))}P with f(x) = a/sqrt(x):")
    print(f"\ta = {params[0]} ± {err[0]}")
    #print(f"\tb = {params[1]} ± {err[1]}")

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
        file_B_field = "[" + filenames_convergence_new[i].split("[")[1].split("]")[0] + "]"
        if len(B_field_names) < 10:
            file_B_field = ""

        B_field_label = ""
        if file_B_field == "[001]":
            B_field_label = hc.B001
        elif file_B_field == "[111]":
            B_field_label = hc.B111
        elif file_B_field == "[011]":
            B_field_label = hc.B011
        else:
            print(f"Warning: Unknown B-field direction {B_field_names}.\n")



        #label = (str(extract_number(filenames_convergence_new[i][7:])) + r"$^3$ lattice sites with " + B_field_label)
        label =  B_field_label

        x = np.asarray(system_amount)
        y = np.asarray(differences[i])

        # plot marker/points on main axes and record the returned Line2D to obtain its color
        main_line, = ax.plot(x, y, "^", label=label)

        # get color from the just-plotted line so the fit uses the same color
        color = main_line.get_color()

        # plot fit curve on main axes using the same color
        fit_line_x = N_fit
        fit_line_y = fit_sqrt(N_fit, *params)
        ax.plot(fit_line_x, fit_line_y, "--", color=color)

       # For the fit: you plotted fit(N_fit, ...) vs 1/N earlier; keep consistent:
        fit_x = 1.0 / np.sqrt(np.asarray(N_fit))
        fit_y = fit_sqrt(N_fit, *params)

        # Plot into inset for overview
        x = 1.0/np.sqrt(x)
        axins.plot(x, y, "^", color=color)
        axins.plot(fit_x, fit_y, "--", color=color)

        # update global y min/max for the inset limits
        all_y_min = min(all_y_min, np.nanmin(y), np.nanmin(fit_y))
        all_y_max = max(all_y_max, np.nanmax(y), np.nanmax(fit_y))

    # tidy up main axes
    ax.set_xlabel( hc.interim_amount_label, fontsize=hc.size_label)
    ax.set_ylabel( hc.deviation_3D_label_interim, fontsize=hc.size_label)
    ax.legend(loc = "center",bbox_to_anchor=(0.85, 0.25), fontsize='small')

    axins.set_xlabel(r"$1/\sqrt{N}$")
    axins.locator_params(axis='x', nbins=7)
    #axins.set_ylabel(r'$\langle \Delta F_x \rangle$ / $\%$')

    # Make the inset show the whole data range (x from main axes and y from aggregated min/max)
    axins.set_xlim(axins.get_xlim())
    # add small padding to avoid data touching the border
    pad = 0.05 * (all_y_max - all_y_min) if (all_y_max > all_y_min) else 0.1
    axins.set_ylim(all_y_min - pad, all_y_max + pad)

    pdfname = pdfname + "_inlet"  # remove "_inverse" suffix for first plot
    outpath = "build/plots/convergence/" + pdfname + ".pdf"
    fig.savefig(outpath, bbox_inches='tight')

    # If you want to replace matched files in the thesis folder
    matching_paths = find_pdfs(pdfname)
    for path in matching_paths:
        print(f"\n\tReplacing plot in master thesis folder: {path}")
        fig.savefig(path, bbox_inches='tight')

    plt.close(fig)


    plt.figure()
    plt.rc('axes', labelsize=size_label)
    #plt.axhline(y = 0, xmin = 0, xmax = 1,linestyle ="dotted", color = "black")
    for i in range(0, len(differences)):

        file_B_field = "[" + filenames_convergence_new[i].split("[")[1].split("]")[0] + "]"
        if(len(B_field_names) < 10):
            file_B_field  = ""
        plt.plot(1/np.sqrt(system_amount ), differences[i ],"^", label=B_field_label)
        #plt.plot(system_amount, differences[i ], color = plt.gca().lines[-1].get_color())
        plt.plot(1/np.sqrt(N_fit), fit_sqrt(N_fit, *params ),"--", color = plt.gca().lines[-1].get_color())#, label=str(extract_number(filenames_convergence_new[i][7:])) + r"$^3$ Particles fit")
    plt.xlabel( "$1/\sqrt{N}$")
    plt.locator_params(axis='x', nbins=7)    #plt.xscale('function', functions=(inverse, inverse))
    plt.ylabel( hc.deviation_3D_label_interim)
    plt.legend()    
    plt.tight_layout()

    pdfname = pdfname[:-8] + "_inverse"
    plt.savefig("build/plots/convergence/" + pdfname + ".pdf")
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



#get rid of pair correlation files



#mask = "pair" not in all_filenames
#all_filenames = all_filenames[mask]
mask = [False] * len(all_filenames)
for i in range(0, len(all_filenames)):
    mask[i] = "pair" not in all_filenames[i] and ".txt" in all_filenames[i] and "rng" not in all_filenames[i] and "no_norm" in all_filenames[i] and "averaged_FID" not in all_filenames[i]
all_filenames = all_filenames[mask]

allready_checked_sys_amount = np.array([])
allready_checked_B_field = np.array([])
allready_checked_sys_size = np.array([])

interim_name = np.array([])
final_name = np.array([])


for file in all_filenames:

    
    system_amount = extract_number(file[:6])
    sys_size = extract_number(file[10:])
    #sort for B_fields
    B_field_name = "[" + file.split("[")[1].split("]")[0] + "]"

    continue_flag = True

    for i in range(0, allready_checked_sys_amount.size):
        if system_amount == allready_checked_sys_amount[i] and B_field_name == allready_checked_B_field[i] and sys_size == allready_checked_sys_size[i]:
            #print(f"Already checked {system_amount} and {B_field_name}")
            continue_flag = False

    if continue_flag == False:
        continue


    mask = [False] * len(all_filenames)
    for i in range(0, len(all_filenames)):
        if(system_amount == extract_number(all_filenames[i][:8]) ) and sys_size == extract_number(all_filenames[i][10:]):
            mask[i] = True
    grouped_filenames = all_filenames[mask]

    #print(f"system amount {system_amount} mask{mask} groups{grouped_filenames}")

    allready_checked_sys_amount = np.append(allready_checked_sys_amount, system_amount)
    allready_checked_B_field = np.append(allready_checked_B_field, B_field_name)
    allready_checked_sys_size = np.append(allready_checked_sys_size, sys_size)

    interim_mask = [False] * len(grouped_filenames)
    for i in range(0, len(grouped_filenames)):
        if("_interim_Correlation_results" in grouped_filenames[i]):
            interim_mask[i] = True

    interim_name = grouped_filenames[interim_mask]

    final_mask = [False] * len(grouped_filenames)
    for i in range(0, len(grouped_filenames)):  
        if("averaged_Correlation" in grouped_filenames[i] and "interim" not in grouped_filenames[i]):
            final_mask[i] = True
            
    final_name = grouped_filenames[final_mask]


    interim_name = np.sort(interim_name)
    final_name = np.sort(final_name)



    mask = [False] * len(interim_name)
    for file in final_name:
        system_size_temp= extract_number(file[10:])
        B_field_temp = "[" + file.split("[")[1].split("]")[0] + "]"
        mask_temp = [False] * len(interim_name)
        for i in range(0, interim_name.size):
            if(system_size_temp == extract_number(interim_name[i][10:]) and B_field_temp in interim_name[i]):
                mask_temp[i] = True
        mask = np.logical_or(mask, mask_temp)

    interim_name = interim_name[mask]
            
    mask = [False] * len(final_name)

    for file in interim_name:
        system_size_temp= extract_number(file[10:])
        B_field_temp = "[" + file.split("[")[1].split("]")[0] + "]"
        mask_temp = [False] * len(final_name)
        for i in range(0, final_name.size):
            if(system_size_temp == extract_number(final_name[i][10:]) and B_field_temp in final_name[i]):
                mask_temp[i] = True
        mask = np.logical_or(mask, mask_temp)
    final_name = final_name[mask]

    if interim_name.size != final_name.size:
        print(f"Warning: Different amount of final and interim files for system amount {system_amount}.\n")
        print(f"Final files: {final_name}\nInterim files: {interim_name}\n")
        continue
    if interim_name.size == 0:
        print(f"No matching files for system amount {system_amount}.\n")
        continue

    plot_correlation_convergence(final_name, interim_name)#, exp_111, t_111)                    plot_interim_difference(file, interim_results_file)
    for i in range(0, interim_name.size):
        plot_interim_difference(final_name[i], interim_name[i])#, exp_111, t_111)                    plot_interim_difference(file, interim_results_file)
    plot_correlation_convergence_sqrt(final_name, interim_name)

