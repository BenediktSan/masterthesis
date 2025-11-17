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
import sys
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

import header_constants as hc


if os.path.exists("build") == False:
    os.mkdir("build")

if os.path.exists("build/plots") == False:
    os.mkdir("build/plots")


if os.path.exists("build/plots/comparison") == False:
    os.mkdir("build/plots/comparison")


#include constants to rescale time

#constants from book 
gamma = 251.662 * 10**6 # rad * µs/T   40.077;//40.069244;//40.077 ;// 1/(T µs)
lattice_const = 2.7315 *10**-10# 5.451;// 2.724;//5.451; //in Angstrom

hbar =const.hbar
mu_0 = const.mu_0

coupling = mu_0 * gamma**2 * hbar**2 / ( 4 * np.pi *lattice_const**3)
factor = hbar**2  * 10**30 
#print(f"Coupling {coupling}\nfactor {factor}\nµ_0 {mu_0}\nhbar {hbar}")


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


def fit(x, a, b):
    return a * x + b

def do_linear_fit(file_name, inv_sys_amount, deviations):

    popt, pcov = curve_fit(fit, inv_sys_amount, deviations)
    perr = np.sqrt(np.diag(pcov))
    print(f"\nFitting deviations from {file_name}")
    print(f"\tFit parameters: \n\ta = ({popt[0] * 1000:.3f} ± {perr[0] * 1000:.3f}) e-3, \n\tb = ({popt[1] * 1000:.3f} ± {perr[1] * 1000:.3f}) e-3")
    print(f"\tInfinite systemm interpolation = {fit(0, *popt) * 1000:.3f} e-3\n")

    return popt, perr

def do_inlet_plot(time, all_correlations, all_err, exp_time, exp_data, filenames_comparison_new, name, special_name, old_name):
    

    #append for nicer figure
    factor = 0.8
    zero = np.zeros( int(len(all_correlations[0])*factor) )

    inlet_correlations = []
    inlet_time = []
    for i in range(len(all_correlations)):
        #inlet_correlations.append(np.append(all_correlations[i], zero))
        #inlet_time.append(np.append(time[0],zero))
        inlet_correlations.append(all_correlations[i])
        inlet_time.append(time[i])
    
    #inlet_correlations.append(all_correlations[0])
    #inlet_correlations.append(all_correlations[-1])
    #inlet_time.append(time[0])
    #inlet_time.append(time[-1])


    size_label = hc.size_label

    fig, ax = plt.subplots()
    plt.rc('axes', labelsize=size_label)

    # horizontal zero line on main axes
    ax.axhline(y=0, xmin=0, xmax=1, linestyle="dotted", color="black")

    # create inset axes (30% x 30% of parent axes, upper right)
    axins = inset_axes(ax, width="65%", height="40%", loc ="upper right")

    # choose log mode:
    use_symlog = False   # set True if you want symmetric log for negative values
    eps = 0#1e-12          # small offset if you prefer shifting data to positive


    axins.set_yscale('log')

    # Plot everything: plot same data into main axes and inset so colors match
    ax.set_prop_cycle(None)    # reset color cycle for main
    axins.set_prop_cycle(None) # reset color cycle for inset so they match
    #axins.set_prop_cycle(None) # reset color cycle for inset so they match
    #axins.plot(exp_time, abs(exp_data), color="black", alpha=0.75)

    for i in range(len(all_correlations)):
        y = all_correlations[i]
        y_for_inset = y.copy()
        y_for_inset = abs(y_for_inset)  # inset is log scale, so take abs value

        if use_symlog:
            # symlog accepts negative / zero values
            pass
        else:
            # log scale needs strictly positive values — either shift or skip non-positive points
            if np.any(y_for_inset <= 0):
                # Option A: shift upward by epsilon (keeps shape but inaccurate for absolute scale)
                y_for_inset = y_for_inset 
                # Option B (commented): skip plotting points <= 0 (creates gaps)
                # mask = y_for_inset > 0
                # axins.plot(time[i][mask], y_for_inset[mask], label=..., ...)
        ax.plot(time[i], all_correlations[i],
                label=str(extract_number(filenames_comparison_new[i][7:])) + r"$^3$ lattice sites " + special_name[i] + old_name[i])
        if(i == 0 or i == len(all_correlations)-1):
            axins.plot(inlet_time[i], abs(inlet_correlations[i]), color = ax.lines[-1].get_color())  # same data (possibly shifted) in inset
    axins.plot(exp_time, abs(exp_data), color=hc.exp_color, alpha=0.75)

    # Plot filled error bands on main axes and inset if desired
    ax.set_prop_cycle(None)
    axins.set_prop_cycle(None)
    for i in range(len(all_correlations)):
        
        y_log = np.asarray(all_correlations[i])
        t = np.asarray(time[i])

        # If all_err is the error in the log-domain (i.e., ± on log|f|), exponentiate bounds:
        upper_lin = np.exp(y_log + all_err[i])
        lower_lin = np.exp(y_log - all_err[i])

        # Fill main (log-domain). You probably have log-domain errors, so fill using log values:
        #ax.fill_between(t, y_log + all_err[i], y_log - all_err[i], alpha=0.3)


        # Fill inset (linear-domain on log scale) using exponentiated bounds
        #axins.fill_between(t, lower_lin, upper_lin, alpha=0.2)

    # experimental data on main axes
    ax.plot(exp_time, exp_data, color="black", alpha=0.75, label="Experimental Data")

    ax.set_xlabel(hc.time_label, fontsize = hc.size_label)
    ax.set_ylabel(hc.FID_label, fontsize = hc.size_label)
    ax.legend(loc = "center",bbox_to_anchor=(0.84, 0.35), fontsize='small')

    #ax.legend(loc='lower left', fontsize='small')

    # optional: remove ticks from inset to reduce clutter
    #axins.set_xticks([])
    axins.set_yticks([1,1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6])

    # optional: draw connectors between inset and zoom region (you can set x1,x2,y1,y2 to zoom)
    # x1, x2 = 6, 9
    # y1, y2 = 1e-3, 1e-1
    # axins.set_xlim(x1, x2)
    # axins.set_ylim(y1, y2)
    # mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")

    pdfname = "Correlation_comparison_" + name + "_inlet"
    fig.savefig("build/plots/comparison/" + pdfname + ".pdf")
    matching_paths = find_pdfs( pdfname)
    for path in matching_paths:
        print(f"\n\tReplacing plot in master thesis folder: {path}")
        plt.savefig(path )
    plt.close()


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



def calc_difference(filenames, comparison_file = ""):

    #load comparison data
    comp_data = np.array([])
    comp_time = np.array([])
    
    #default to experimental data
    if comparison_file == "":

        #choose exp data
        if "[001]" in filenames[0]:
            comp_data = exp_100
            comp_time = t_100
        elif "[011]" in filenames[0]:
            comp_data = exp_110
            comp_time = t_110        
        elif "[111]" in filenames[0]:
            comp_data = exp_111
            comp_time = t_111
    else:
        comp_data, _, comp_time = np.genfromtxt("build/daten/"+comparison_file,skip_header= 1, unpack = True)

    #get data
    all_correlations, all_err, t = np.genfromtxt("build/daten/"+filenames[0],skip_header= 1, unpack = True)
    all_correlations = np.array([all_correlations])
    all_err = np.array([all_err])

    #conversion of t for interpolation
    t_tilde = t /(hbar )*10**-30 * 10**6 

        
    #And now we also do a retrospective correction to the lattice constant, because I used a wrong one in th e simulation
    new_a = 2.72325
    my_a = 2.7315
    factor_a = (new_a/my_a)**(-3)
    #print( f"factor a = {factor_a}")
    t_tilde = t_tilde / factor_a
    
    #normlaise 0th value in array
    all_correlations = all_correlations / all_correlations[0][0]

    amount_of_files = len(filenames)

    for i in range(1, amount_of_files):

        correlation, err, t_new = np.genfromtxt("build/daten/"+filenames[i],skip_header= 1, unpack = True)
        if(len(t) != len(t_new)):
            print(f"##########ERROR SAME AMOUNT OF INITIAL CONDITIONS BUT DIFFERENT STEP AMOUNT########")
            continue
        
        #normalize in case something is not normalized
        #err sucks then tho
        correlation = correlation / correlation[0]

        all_correlations = np.append(all_correlations, [correlation], axis = 0)
        all_err = np.append(all_err, [err], axis = 0)


    #start calculation

    interpolated_values = np.zeros(len(t_tilde))
    for i in range(0, len(t_tilde)):
        interpolated_values[i] = interpolate(t_tilde[i], comp_time, comp_data)

    #plot_interpolation(comp_time, comp_data, t_tilde, interpolated_values)


    deviation = np.zeros(amount_of_files)
    std_dev = np.zeros(amount_of_files)

    for i in range(0, amount_of_files):
        delta = 0
        step_amount = len(t_tilde)
        I = sum(abs(comp_data))
        for j in range(0,step_amount):
            delta += abs((interpolated_values[j] - all_correlations[i][j])/ I)
        #maybe even take the square root
        deviation[i] =  delta / step_amount

    for i in range(0, amount_of_files):
        delta_squared = 0
        step_amount = len(t_tilde)
        I = sum(abs(comp_data))
        for j in range(0,step_amount):
            delta_squared += (abs((all_correlations[i][j] - interpolated_values[j])/ I) - deviation[i])**2
        #maybe even take the square root
        std_dev[i] =np.sqrt( delta_squared / (step_amount-1))


    print(filenames[0])
    print(deviation[0])


    return deviation, std_dev


def plot_interpolation(t_exp, exp_data, t, interp_val):
    size_label = hc.size_label
    plt.figure()
    plt.rc('axes', labelsize=size_label)
    plt.axhline(y = 0, xmin = 0, xmax = 1,linestyle ="dotted", color = "black")
    plt.plot(t_exp, exp_data, label=r"exp Data")
    plt.plot(t, interp_val, label=r"interpolated Values")
    #plt.plot(t, M_x_comparison[1], label=r"M_x " + filenames_comparison[1][-31:-29])
    #plt.plot(t, M_x_comparison[2], label=r"M_x " + filenames_comparison[2][-32:-29])
    plt.xlabel( hc.time_label )
    plt.ylabel( hc.FID_label )
    plt.legend()
    plt.tight_layout()
    plt.savefig("build/plots/Interpolation.pdf")
    pdfname = "Interpolation"
    plt.savefig("build/plots/" + pdfname + ".pdf")
    #check if plot is in master thesis folder and replace it if it is
    matching_paths = find_pdfs( pdfname)
    for path in matching_paths:
        print(f"\n\tReplacing plot in master thesis folder: {path}")
        plt.savefig(path )
    plt.close()

def plot_comparison(filenames_comparison):



    if(len(filenames_comparison) < 2):
        return
    print(f"\nstarted plotting for \n{filenames_comparison}\n\n")


    old_name = []

    for file in filenames_comparison:
        if "_old" in file:
            old_name.append("old")
        else:
            old_name.append("")

    old_name = np.asarray(old_name)



    #get data
    all_correlations, all_err, t = np.genfromtxt("build/daten/"+filenames_comparison[0],skip_header= 1, unpack = True)
    all_correlations = np.array([all_correlations])
    all_err = np.array([all_err])


    #normlaise 0th value in array
    all_correlations = all_correlations / all_correlations[0][0]

    for i in range(1,len(filenames_comparison)):

        correlation, err, t_new = np.genfromtxt("build/daten/"+filenames_comparison[i],skip_header= 1, unpack = True)
        if(len(t) != len(t_new)):
            print(f"##########ERROR SAME AMOUNT OF INITIAL CONDITIONS BUT DIFFERENT STEP AMOUNT########")
            continue
        #normalize in case something is not normalized
        #err sucks then tho
        correlation = correlation / correlation[0]

        all_correlations = np.append(all_correlations, [correlation], axis = 0)
        all_err = np.append(all_err, [err], axis = 0)


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
    #print(f"\t\t450 steps are {450 * t_tilde[1]} {t_tilde[-1]} µs for 100µs {200/t_tilde[1]} steps are needed")  

    time = []

    for i in range(0, len(all_correlations)):
        
        time.append(t_tilde)

    time = np.asarray(time)

    for i in range(0, len(time)):
        if("_old" in filenames_comparison[i]):
            time[i] /= 2



    #name for file
    filenames_comparison_new = np.array([])
    filenames_comparison_new = np.append(filenames_comparison_new,filenames_comparison)


    #first removing lido3 if in file
    if("lido3" in str(filenames_comparison_new)):
        for i in range(0,filenames_comparison_new.size):
            filenames_comparison_new[i] = filenames_comparison_new[i].replace("lido3_", "")



    name = str(extract_number(filenames_comparison_new[0][:7])) +"sys"
    for file in filenames_comparison_new:
        name += "_" + str(extract_number(file[7:])) + "P"

    #implies only one kind of magentic field
    B_field_name = "[" + filenames_comparison_new[0].split("[")[1].split("]")[0] + "]"
    name = B_field_name + "_" + name

    #get experimental data
    exp_data = np.array([])
    exp_time = np.array([])

    #choose exp data
    if "[001]" in filenames_comparison[0]:
        exp_data = exp_100
        exp_time = t_100
    elif "[011]" in filenames_comparison[0]:
        exp_data = exp_110
        exp_time = t_110        
    elif "[111]" in filenames_comparison[0]:
        exp_data = exp_111
        exp_time = t_111


    special_name = []
    for file in filenames_comparison_new:
        spc_name = ""
        if "rng" in file or "no_norm" not in file:
            spc_name += " with "
        if "rng" in file:
            spc_name += "Gaussian i.c. "
        if "no_norm" not in file:
            spc_name += r"N_2" + "used"
        special_name.append( spc_name)


    size_label = hc.size_label
    plt.figure()
    plt.rc('axes', labelsize=size_label)
    plt.axhline(y = 0, xmin = 0, xmax = 1,linestyle ="dotted", color = "black")
    plt.gca().set_prop_cycle(None)
    for i in range(0, len(all_correlations)):
        plt.plot(time[i], all_correlations[i], label= str(extract_number(filenames_comparison_new[i][7:])) + r"$^3$ lattice sites " + special_name[i] + old_name[i])
    plt.gca().set_prop_cycle(None)
    for i in range(0, len(all_correlations)):
        plt.fill_between(time[i], all_correlations[i] + all_err[i], all_correlations[i] - all_err[i], alpha = 0.3)
    plt.plot(exp_time, exp_data, color = hc.exp_color, alpha = 0.75, label = "Experimental data")
    #plt.plot(t, M_x_comparison[0], label=r"M_x " + filenames_comparison[0][-31:-29])
    #plt.plot(t, M_x_comparison[1], label=r"M_x " + filenames_comparison[1][-31:-29])
    #plt.plot(t, M_x_comparison[2], label=r"M_x " + filenames_comparison[2][-32:-29])
    plt.xlabel( hc.time_label )
    plt.ylabel( hc.FID_label )
    plt.legend()
    plt.tight_layout()
    pdfname = "Correlation_comparison_" + name
    plt.savefig("build/plots/comparison/" + pdfname + ".pdf")
    #check if plot is in master thesis folder and replace it if it is
    matching_paths = find_pdfs( pdfname)
    for path in matching_paths:
        print(f"\n\tReplacing plot in master thesis folder: {path}")
        plt.savefig(path )
    plt.close()

    do_inlet_plot(time, all_correlations, all_err, exp_time, exp_data, filenames_comparison_new, name, special_name, old_name)


    #add additional log plot
    plt.figure()
    plt.rc('axes', labelsize=size_label)
    plt.axhline(y = 0, xmin = 0, xmax = 1,linestyle ="dotted", color = "black")
    plt.gca().set_prop_cycle(None)
    for i in range(0, len(all_correlations)):
        plt.plot(time[i], abs(all_correlations[i]), label= str(extract_number(filenames_comparison_new[i][7:])) + r"$^3$ lattice sites " + special_name[i] + old_name[i])
    plt.gca().set_prop_cycle(None)
    #for i in range(0, len(all_correlations)):
        #plt.fill_between(time[i], all_correlations[i] + all_err[i], all_correlations[i] - all_err[i], alpha = 0.3)
    plt.plot(exp_time, abs(exp_data), color = hc.exp_color, alpha = 0.75, label = "Experimental data")
    #plt.plot(t, M_x_comparison[0], label=r"M_x " + filenames_comparison[0][-31:-29])
    #plt.plot(t, M_x_comparison[1], label=r"M_x " + filenames_comparison[1][-31:-29])
    #plt.plot(t, M_x_comparison[2], label=r"M_x " + filenames_comparison[2][-32:-29])
    plt.xlabel( hc.time_label )
    plt.ylabel( r"$|$" + hc.FID_label + r"$|$" )
    plt.yscale("log")
    plt.tight_layout()
    #plt.ylim([10**(-6),None])    
    plt.legend()
    pdfname = "Correlation_comparison_" + name + "_log"
    plt.savefig("build/plots/comparison/" + pdfname + ".pdf")
    #check if plot is in master thesis folder and replace it if it is
    matching_paths = find_pdfs( pdfname)
    for path in matching_paths:
        print(f"\n\tReplacing plot in master thesis folder: {path}")
        plt.savefig(path )
    plt.close()


def plot_comparison_difference(filenames_comparison):
    
    if(len(filenames_comparison) < 2):
        return
    #implies only one kind of magentic field
    B_field_name = "[" + filenames_comparison[0].split("[")[1].split("]")[0] + "]"

    differences, std_dev = calc_difference(filenames_comparison)

    #extract lido3 to get particle amount
    #name for file
    filenames_comparison_new = np.array([])
    filenames_comparison_new = np.append(filenames_comparison_new,filenames_comparison)


    if("lido3" in str(filenames_comparison_new)):
        for i in range(0,filenames_comparison_new.size):
            filenames_comparison_new[i] = filenames_comparison_new[i].replace("lido3_", "")


    particle_amount = np.array([extract_number(filenames_comparison_new[0][7:])])
    for i in range(1, len(filenames_comparison_new)):
        particle_amount= np.append(particle_amount,[extract_number(filenames_comparison_new[i][7:])])
    #print(particle_amount)
    particle_amount_cubed = np.power(particle_amount,3)


    special_name = []
    for file in filenames_comparison_new:
        spc_name = ""
        if "rng" in file or "no_norm" not in file:
            spc_name += " with "
        if "rng" in file:
            spc_name += "Gaussian i.c. "
        if "no_norm" not in file:
            spc_name += r"N_2" + "used"
        special_name.append( spc_name)

    ticks, tick_labels = hc.gain_ticks_inverse_system_size(particle_amount)

    #do fit for deviations
    params, err = do_linear_fit(filenames_comparison[0], 1/particle_amount_cubed, differences)


    size_label = hc.size_label
    plt.figure()
    plt.rc('axes', labelsize=size_label)
    #plt.axhline(y = 0, xmin = 0, xmax = 1,linestyle ="dotted", color = "black")
    for i in range(0, len(particle_amount)):

        plt.plot(1/particle_amount_cubed[i], differences[i],"^", label= str(particle_amount[i]) + r"$^3$ lattice sites" + special_name[i])
        plt.errorbar(1/particle_amount_cubed[i], differences[i], yerr = std_dev[i], fmt ="^", color = plt.gca().lines[-1].get_color(), capsize=3)

    if(B_field_name != "[001]"):
        plt.plot(1/particle_amount_cubed, fit(1/particle_amount_cubed, *params), linestyle="dashed", color = hc.fit_color, label=f"Linear fit")

    #plt.plot(t, M_x_comparison[0], label=r"M_x " + filenames_comparison[0][-31:-29])
    #plt.plot(t, M_x_comparison[1], label=r"M_x " + filenames_comparison[1][-31:-29])
    #plt.plot(t, M_x_comparison[2], label=r"M_x " + filenames_comparison[2][-32:-29])
    plt.xlabel( hc.inverse_system_size_label )
    plt.ylabel( hc.deviation_3D_label )
    plt.xticks(ticks, tick_labels)
    plt.legend(loc='upper right', fontsize = "large", bbox_to_anchor=(0.94, 1))
    plt.tight_layout()
    pdfname = "Correlation_deviation_" + B_field_name + "_" + str(extract_number(filenames_comparison[0][:7])) + "sys"
    plt.savefig("build/plots/comparison/" + pdfname + ".pdf")
    #check if plot is in master thesis folder and replace it if it is
    matching_paths = find_pdfs( pdfname)
    for path in matching_paths:
        print(f"\n\tReplacing plot in master thesis folder: {path}")
        plt.savefig(path )    
    plt.close()
    
def plot_comparison_difference_maximum_system_size(filenames_comparison):
    
    if(len(filenames_comparison) < 3):
        return
    
    filenames_comparison.sort()

    #intorduce file used for compariosn
    comparison_file = filenames_comparison[-1]
    filenames_comparison = np.delete(filenames_comparison,-1,0)

    #implies only one kind of magentic field
    B_field_name = "[" + filenames_comparison[0].split("[")[1].split("]")[0] + "]"

    differences, std_dev = calc_difference(filenames_comparison, comparison_file)

    #extract lido3 to get particle amount
    #name for file
    filenames_comparison_new = np.array([])
    filenames_comparison_new = np.append(filenames_comparison_new,filenames_comparison)


    if("lido3" in str(filenames_comparison_new)):
        for i in range(0,filenames_comparison_new.size):
            filenames_comparison_new[i] = filenames_comparison_new[i].replace("lido3_", "")


    particle_amount = np.array([extract_number(filenames_comparison_new[0][7:])])
    for i in range(1, len(filenames_comparison_new)):
        particle_amount= np.append(particle_amount,[extract_number(filenames_comparison_new[i][7:])])
    #print(particle_amount)
    particle_amount_cubed = np.power(particle_amount,3)


    special_name = []
    for file in filenames_comparison_new:
        spc_name = ""
        if "rng" in file or "no_norm" not in file:
            spc_name += " with "
        if "rng" in file:
            spc_name += "Gaussian i.c. "
        if "no_norm" not in file:
            spc_name += r"N_2" + "used"
        special_name.append( spc_name)

    ticks, tick_labels = hc.gain_ticks_inverse_system_size(particle_amount)

    max_sys_name = str(extract_number(grouped_filenames[-1][10:]))
    max_sys = extract_number(grouped_filenames[-1][10:])

    print(differences)
    print("\nFitting deviations compared to maximum system size:")
    params, err = do_linear_fit(filenames_comparison[0], 1/particle_amount_cubed, differences)

    size_label = hc.size_label
    plt.figure()
    #plt.title(r"Deviation in regard to the FID for $L=$" + max_sys_name )
    plt.rc('axes', labelsize=size_label)
    #plt.axhline(y = 0, xmin = 0, xmax = 1,linestyle ="dotted", color = "black")
    for i in range(0, len(particle_amount)):

        plt.plot(1/particle_amount_cubed[i], differences[i],"^", label= str(particle_amount[i]) + r"$^3$ lattice sites" + special_name[i])
        #plt.errorbar(1/particle_amount_cubed[i], differences[i], yerr = std_dev[i], fmt ="^", color = plt.gca().lines[-1].get_color(), capsize=3)
    if(B_field_name != "[001]"):
        plt.plot(1/particle_amount_cubed, fit(1/particle_amount_cubed, *params), linestyle="dashed", color = hc.fit_color, label=f"Linear fit")

    #plt.plot(t, M_x_comparison[0], label=r"M_x " + filenames_comparison[0][-31:-29])
    #plt.plot(t, M_x_comparison[1], label=r"M_x " + filenames_comparison[1][-31:-29])
    #plt.plot(t, M_x_comparison[2], label=r"M_x " + filenames_comparison[2][-32:-29])
    plt.xlabel( hc.inverse_system_size_label)
    plt.ylabel( hc.deviation_3D_label )
    plt.xticks(ticks, tick_labels)
    plt.legend(loc = "best", fontsize = "large")
    plt.tight_layout()
    pdfname = "Correlation_deviation_" + B_field_name + "_" + str(extract_number(filenames_comparison[0][:7])) + "sys_regard_to_Maximum_" + max_sys_name + "P"
    plt.savefig("build/plots/comparison/" + pdfname + ".pdf")
    #check if plot is in master thesis folder and replace it if it is
    matching_paths = find_pdfs( pdfname)
    for path in matching_paths:
        print(f"\n\tReplacing plot in master thesis folder: {path}")
        plt.savefig(path )    
    plt.close()
    

'''
# dependencies for a comparison plot
do_plot = False
t_comparison = np.array([0])
C_x_comparison = np.array([0])
#M_y_comparison = np.array([0])
#M_z_comparison = np.array([0])
filenames_comparison = np.array([0])
'''


all_filenames = np.array([])

for file in os.listdir("build/daten"):
    all_filenames = np.append(all_filenames,file)

#get rid of pair correlation files



#mask = "pair" not in all_filenames
#all_filenames = all_filenames[mask]
mask = [False] * len(all_filenames)
for i in range(0, len(all_filenames)):
    mask[i] = "pair" not in all_filenames[i] and "interim" not in all_filenames[i] and "variance" not in all_filenames[i] and ".txt" in all_filenames[i] 
all_filenames = all_filenames[mask]

# just so its not emoty and i can append
allready_checked_sys_amount = np.array([])
allready_checked_B_field = np.array([])


for file in all_filenames:
    
    system_amount = extract_number(file[:6])
    #sort for B_fields
    B_field_name = "[" + file.split("[")[1].split("]")[0] + "]"

    continue_flag = True

    for i in range(0, allready_checked_sys_amount.size):
        if system_amount == allready_checked_sys_amount[i] and B_field_name == allready_checked_B_field[i]:
            #print(f"Already checked {system_amount} and {B_field_name}")
            continue_flag = False

    if continue_flag == False:
        continue

    mask = [False] * len(all_filenames)
    for i in range(0, len(all_filenames)):
        if(system_amount == extract_number(all_filenames[i][:8]) and B_field_name in all_filenames[i]):
            mask[i] = True
    grouped_filenames = all_filenames[mask]

    #print(f"system amount {system_amount} mask{mask} groups{grouped_filenames}")

    allready_checked_sys_amount = np.append(allready_checked_sys_amount, system_amount)
    allready_checked_B_field = np.append(allready_checked_B_field, B_field_name)

    #sort them foe aesthetic and for comparioson to maximum size
    grouped_filenames.sort()

    plot_comparison(grouped_filenames)
    plot_comparison_difference(grouped_filenames)
    plot_comparison_difference_maximum_system_size(grouped_filenames)



'''

        
if(do_plot == True):
    plot_comparison(filenames_comparison)
        
print("Dateianzahl: ", file_amount)
'''
