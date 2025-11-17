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
import copy
import sys

import header_constants as hc

if os.path.exists("build") == False:
    os.mkdir("build")

if os.path.exists("build/plots") == False:
    os.mkdir("build/plots")


if os.path.exists("build/plots/hybrid_approach") == False:
    os.mkdir("build/plots/hybrid_approach")


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


def fit(t, C, omega, psi, gamma):
    return C * np.exp(-gamma*t) * np.cos(omega*t + psi)

def calc_fit_params(name,t, C):
    #initial guess
    C_0 = 0.6
    omega_0 = 160 #MHz
    psi_0 = -0.5
    gamma_0 = 0.052 #1/µs

    popt, pcov = curve_fit(fit, t, C, p0=[C_0, omega_0, psi_0, gamma_0])

    print(f"\n\tFit parameters for {name}:\n\tC = {popt[0]}\n\tomega = {popt[1]} MHz\n\tpsi = {popt[2]} rad\n\tgamma = {popt[3]} MHz\n")

    return popt, pcov





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
                timo_correlations = timo_data['CaF-' + new_B_field_name +' p.-c. contr. G^xx from nl-spinDMFT'][1:]
                break
    return timo_t, timo_correlations

def calc_hybrid_FID( timo_correlations,timo_t, pair_and_auto_correlations, t_pair):
    """
    Calculate the hybrid FID from Timo's data and the pair correlations.
    
    :param timo_t: Time values from Timo's data
    :param timo_correlations: Correlation values from Timo's data
    :param pair_and_auto_correlations: Pair and auto correlations from the current data
    :return: Hybrid FID
    """



    # cut array so both have approxiametely the same t_max
    if(timo_t[-1] >= t_pair[-1]):

        timo_t = timo_t[timo_t <= t_pair[-1]]
        timo_correlations = timo_correlations[:len(timo_t)]
    else:

        t_pair = t_pair[t_pair <= timo_t[-1]]
        pair_and_auto_correlations = pair_and_auto_correlations[:, :len(t_pair)]
       
        

    #allready add up timos correlations to gain hybrid FID
    hybrid_FID = np.zeros(len(timo_t))
    for i in range(0, len(timo_correlations)):
        hybrid_FID += timo_correlations[i]

    #interpolate values 
    #We use Gräßers data as base for the interpolation. 
    #this is because my data has finder time steps, therefore the interpolation is more accurate
    #On the other hand thhis means that we have fewer values to calculate the deviation

    hybrid_time = copy.deepcopy(timo_t)
    to_be_interpolated = copy.deepcopy(pair_and_auto_correlations)

    if(len(timo_t) > len(t_pair)):
        print(f"\n\t\tTimo's data has more time steps than the pair correlations. Interpolating Timo's data to the pair correlations.")


    for i in range(0, len(to_be_interpolated)):
        temp_corr = np.zeros(len(hybrid_time))
        for j in range(0, len(hybrid_time) ):
            temp_corr[j] = interpolate(hybrid_time[j], t_pair, to_be_interpolated[i])
        hybrid_FID += temp_corr


    ##normalize FID
    #norm = hybrid_FID[0]
    #hybrid_FID = hybrid_FID / norm


    return hybrid_FID, hybrid_time



def calc_deviation(exp_data, exp_time, timo_correlations, timo_t, pair_and_auto_correlations, t_pair):
    """
    Calculate the deviation between Timo's data and the pair correlations.
    
    :param timo_t: Time values from Timo's data
    :param timo_correlations: Correlation values from Timo's data
    :param pair_and_auto_correlations: Pair and auto correlations from the current data
    :return: Deviation and standard deviation
    """

    ##add up correlations to gain hybrid FID

    hybrid_FID, hybrid_time = calc_hybrid_FID( timo_correlations, timo_t, pair_and_auto_correlations, t_pair)



    #interpolate experimental data to hybrid time
    interpolated_exp_values = np.zeros(len(hybrid_time))
    for i in range(0, len(hybrid_time)):
        interpolated_exp_values[i] = interpolate(hybrid_time[i], exp_time, exp_data)


    deviation = 0
    std_dev = 0

    step_amount = len(hybrid_time)

    delta = 0
    delta_squared = 0

    for j in range(0,step_amount):
        I = sum(abs(interpolated_exp_values))
        delta += abs(interpolated_exp_values[j] - hybrid_FID[j]) / I * 100
        delta_squared += abs((hybrid_FID[j] - interpolated_exp_values[j])/ I * 100)**2
    #maybe even take the square root
    deviation =  delta / step_amount
    std_dev = np.sqrt( delta_squared / (step_amount-1))


    return deviation, std_dev



def calc_different_cut_offs(exp_data, exp_time, pair_filename):



    B_field_name =  "[" + pair_filename.split("[")[1].split("]")[0] + "]"



    data = read_txt_file("build/daten/"+pair_filename)  

    
    B_field_name =  "[" + pair_filename.split("[")[1].split("]")[0] + "]"

    #find Gräßers data



    timo_t, timo_correlations = find_data_for_comparison(pair_filename)
    #read in my data

    pair_correlations = np.array([])


    if(len(data) == 0):
        return  
    
    #convert to np.array
    pair_and_auto_correlations = np.array([data[0]], dtype=np.float32)
    for i in range(1, len(data)):
        pair_and_auto_correlations = np.append(pair_and_auto_correlations,np.array([data[i]],dtype=np.float32), axis = 0)


    #extract t
    t_pair = pair_and_auto_correlations[-1]
    pair_and_auto_correlations = np.delete(pair_and_auto_correlations,-1,0)

    #calculate sum of pair correlations
    sum_of_all_correlations = pair_and_auto_correlations[0]
    for i in range(1, len(pair_correlations)):
        sum_of_all_correlations = sum_of_all_correlations + pair_and_auto_correlations[i]


    #if B_field_name == "[011]":
    #    print(f"\t\tReordering pair correlations for B field {B_field_name}")
    #    temp1 = pair_and_auto_correlations[6].copy()
    #    temp2 = pair_and_auto_correlations[7].copy()
    #    temp3 = pair_and_auto_correlations[8].copy()
    #    temp4 = pair_and_auto_correlations[9].copy()
#
    #    pair_and_auto_correlations[6] = temp2
    #    pair_and_auto_correlations[9] = temp1
    #    pair_and_auto_correlations[7] = temp3
    #    pair_and_auto_correlations[8] = temp4

    if B_field_name == "[011]":
        print(f"\t\tReordering pair correlations for B field {B_field_name}")

        pair_and_auto_correlations[6] += pair_and_auto_correlations[7]
        pair_and_auto_correlations = np.delete( pair_and_auto_correlations,7,0)


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
    #print(f"\n\t\tstepsize for 

    if "old" in pair_filename:
        t_tilde = t_tilde / 2
    t_pair = t_tilde

    ##normalizing
    
    norm =  sum_of_all_correlations[0]

    #norm = 5 ** 3 * 8
    sum_of_all_correlations = sum_of_all_correlations / norm
    pair_and_auto_correlations = pair_and_auto_correlations / norm


    ###start introducing a cut off to Timos QM data and calculate all the deviations


    QM_correlation_amount = len(timo_correlations)
    
    deviation = np.zeros(QM_correlation_amount)
    std_dev = np.zeros(QM_correlation_amount)

    for cut_off in range(0, QM_correlation_amount):
        if(cut_off >= len(pair_and_auto_correlations)):
            print(f"\n\t\tCut off {cut_off} is larger than the amount of pair correlations {len(pair_and_auto_correlations)}")
            deviation[cut_off] = float("nan")
            std_dev[cut_off] = float("nan")
            continue
        
        ##CHECK IF THIS IS THE CORRECT INDEXINGfloat("nan")
        deviation_new, std_dev_new = calc_deviation(exp_data, exp_time,  timo_correlations[:cut_off], timo_t, pair_and_auto_correlations[cut_off:], t_pair)
        deviation[cut_off] = deviation_new
        std_dev[cut_off] = std_dev_new


    return deviation, std_dev




def plot_hybrid_approach(pair_filename):

    print(f"\nstarted plotting for \n{pair_filename}\n\n")
    ###Choose exp data

    B_field_name =  "[" + pair_filename.split("[")[1].split("]")[0] + "]"



    exp_data = np.array([])
    exp_time = np.array([])

    #choose exp data
    if "[001]" in pair_filename:
        exp_data = exp_100
        exp_time = t_100
    elif "[011]" in pair_filename:
        exp_data = exp_110
        exp_time = t_110        
    elif "[111]" in pair_filename:
        exp_data = exp_111
        exp_time = t_111

    #read in my data
    pair_correlations = np.array([])
    data = read_txt_file("build/daten/"+pair_filename)
    if(len(data) == 0):
        print(f"\n\t\tNo data found for {pair_filename}\n")
        return
    #convert to np.array
    pair_and_auto_correlations = np.array([data[0]], dtype=np.float32)
    for i in range(1, len(data)):
        pair_and_auto_correlations = np.append(pair_and_auto_correlations, np.array([data[i]], dtype=np.float32), axis=0)
    #extract t
    t_pair = pair_and_auto_correlations[-1]
    pair_and_auto_correlations = np.delete(pair_and_auto_correlations, -1, 0)

    #calculate sum of pair correlations
    sum_of_all_correlations = pair_and_auto_correlations[0]
    for i in range(1, len(pair_correlations)):
        sum_of_all_correlations = sum_of_all_correlations + pair_and_auto_correlations[i]


    norm = sum_of_all_correlations[0]

    sum_of_all_correlations = sum_of_all_correlations / norm
    pair_and_auto_correlations = pair_and_auto_correlations / norm

    if B_field_name == "[011]":
        print(f"\t\tReordering pair correlations for B field {B_field_name}")
        temp1 = pair_and_auto_correlations[6].copy()
        temp2 = pair_and_auto_correlations[7].copy()
        temp3 = pair_and_auto_correlations[8].copy()
        temp4 = pair_and_auto_correlations[9].copy()

        pair_and_auto_correlations[6] = temp2
        pair_and_auto_correlations[9] = temp1
        pair_and_auto_correlations[7] = temp3
        pair_and_auto_correlations[8] = temp4




    #find Gräßers data for comparison
    t_timo, timo_correlations = find_data_for_comparison(pair_filename)


    #calc all the deviations to identify the optimal cut_off
    deviation, std_dev = calc_different_cut_offs(exp_data, exp_time, pair_filename)


    #very arbtitary stuff happening here
    #factor 1/2 means i miss a factor 2 in my coupling???????????
    t_pair = t_pair /(hbar )*10**-30 * 10**6 

        
    #And now we also do a retrospective correction to the lattice constant, because I used a wrong one in th e simulation
    new_a = 2.72325
    my_a = 2.7315
    factor_a = (new_a/my_a)**(-3)
    #print( f"factor a = {factor_a}")
    t_pair = t_pair / factor_a
    
    #print(f"\t\tmin von exp data t = {t_exp[FID_exp.argmin()]}\n\t\tmin von t_tilde t = {t_tilde[C.argmin()]}\n\t\tmin von t_tilde_2 t = {t_tilde_2[C.argmin()]}")
    #print(f"\n\t\tstepsize for 

    if "old" in pair_filename:
        t_pair = t_pair / 2


    optimal_cut_off = np.where(deviation == min(deviation))
    if(len(optimal_cut_off)>1):
        print(f"\n\t\tMultiple optimal cut offs found for {pair_filename} at {optimal_cut_off}\n")
    optimal_cut_off = int(optimal_cut_off[0][0])

    print(deviation)
    print(f"\n\t\tNumber of QM-Correlations: {len(timo_correlations)} and Number of classical correlations: {len(pair_and_auto_correlations)}")
    print(f"\n\t\tOptimal cut off is {optimal_cut_off} with a deviation of {deviation[optimal_cut_off]} and a std dev of {std_dev[optimal_cut_off]}\n")

    #calculate otimal FID and withc slitghtly suboptimal FIDs

    optimal_FID, hybrid_t = calc_hybrid_FID(timo_correlations[:optimal_cut_off], t_timo, pair_and_auto_correlations[optimal_cut_off:], t_pair)


    suboptimal_factor1 = 0
    suboptimal_factor2 = 0

    if(optimal_cut_off == 0):
        suboptimal_factor1 = 1
        suboptimal_factor2 = 2
    elif(optimal_cut_off == len(timo_correlations) - 1):
        suboptimal_factor1 = -1
        suboptimal_factor2 = -2
    else:
        suboptimal_factor1 = 1
        suboptimal_factor2 = -1    
    

    #do a fit for hybrid FID
    params, std_dev = calc_fit_params(pair_filename, hybrid_t, optimal_FID)



    suboptimal_FID_1, _ = calc_hybrid_FID(timo_correlations[:optimal_cut_off+suboptimal_factor1], t_timo, pair_and_auto_correlations[optimal_cut_off+suboptimal_factor1:], t_pair)
    suboptimal_FID_2, _ = calc_hybrid_FID(timo_correlations[:optimal_cut_off+suboptimal_factor2], t_timo, pair_and_auto_correlations[optimal_cut_off+suboptimal_factor2:], t_pair)


    #get complete QM and classical data

    QM_FID, _ = calc_hybrid_FID(timo_correlations, t_timo, np.array([np.zeros(len(t_pair))]), t_pair)
    classical_FID, _ = calc_hybrid_FID(np.zeros(len(t_timo)), t_timo, pair_and_auto_correlations, t_pair)



    size_label = hc.size_label

    plt.figure()
    plt.rc('axes', labelsize=size_label)
    plt.axhline(y = 0, xmin = 0, xmax = 1,linestyle ="dotted", color = "black")
    #converting into µs and adding all factors neglected in simulation
    plt.plot(hybrid_t , optimal_FID,color = hc.optimal_FID, label=r"Optimal hybrid FID with $c = " + f"{optimal_cut_off}" + r"$")
    #plt.plot(hybrid_t , suboptimal_FID_1,color = "red", label=r"Suboptimal $G_c^{xx}$ with $c = " + f"{optimal_cut_off + suboptimal_factor1}" + r"$")
    #plt.plot(hybrid_t , suboptimal_FID_2,color = "green" ,label=r"Suboptimal $G_c^{xx}$ with $c = " + f"{optimal_cut_off + suboptimal_factor2}" + r"$")
    plt.plot(hybrid_t, QM_FID, color = hc.timo_color_FID, label=r"Gräßer et al.")
    plt.plot(hybrid_t, classical_FID, color = hc.FID_color, label=r"Semi-classical FID")
    plt.plot(hybrid_t , optimal_FID,color = hc.optimal_FID)
    plt.plot(exp_time, exp_data, color = hc.exp_color, label=r"Experimental data")
    #plt.plot( hybrid_t, fit(hybrid_t, *params), linestyle ="dashed", color = hc.fit_color, label = r"Fit: $C e^{-\gamma t} \cos(\omega t + \psi)$")
    plt.xlabel( hc.time_label)
    plt.ylabel( hc.FID_label)
    plt.legend()
    plt.tight_layout()
    pdfname = "hybrid_" + pair_filename[:-4] 
    plt.savefig("build/plots/hybrid_approach/" + pdfname + ".pdf")
    #check if plot is in master thesis folder and replace it if it is
    matching_paths = find_pdfs( pdfname)
    for path in matching_paths:
        print(f"\n\tReplacing plot in master thesis folder: {path}")
        plt.savefig(path )
    plt.close()

    plt.figure()
    plt.rc('axes', labelsize=size_label)
    plt.axhline(y = 0, xmin = 0, xmax = 1,linestyle ="dotted", color = "black")
    #converting into µs and adding all factors neglected in simulation
    plt.plot(hybrid_t , abs(optimal_FID),color = hc.optimal_FID, label=r"Optimal hybrid FID with $c = " + f"{optimal_cut_off}" + r"$")
    #plt.plot(hybrid_t , suboptimal_FID_1,color = "red", label=r"Suboptimal $G_c^{xx}$ with $c = " + f"{optimal_cut_off + suboptimal_factor1}" + r"$")
    #plt.plot(hybrid_t , suboptimal_FID_2,color = "green" ,label=r"Suboptimal $G_c^{xx}$ with $c = " + f"{optimal_cut_off + suboptimal_factor2}" + r"$")
    plt.plot(hybrid_t, abs(QM_FID), color = hc.timo_color_FID, label=r"Gräßer et al.")
    plt.plot(hybrid_t, abs(classical_FID), color = hc.FID_color, label=r"Semi-classical FID")
    plt.plot(hybrid_t , abs(optimal_FID),color = hc.optimal_FID)
    plt.plot(exp_time, abs(exp_data), color = hc.exp_color, label=r"Experimental data")
    plt.yscale("log")
    plt.xlabel( r'$t$ / $µs$')
    plt.ylabel( r'$FID$ ' )
    plt.legend()
    plt.tight_layout()
    pdfname = "hybrid_" + pair_filename[:-4] + "_log"
    plt.savefig("build/plots/hybrid_approach/" + pdfname + ".pdf")
    #check if plot is in master thesis folder and replace it if it is
    matching_paths = find_pdfs( pdfname)
    for path in matching_paths:
        print(f"\n\tReplacing plot in master thesis folder: {path}")
        plt.savefig(path )
    plt.close()

    plt.figure()
    plt.rc('axes', labelsize=size_label)
    plt.axhline(y = 0, xmin = 0, xmax = 1,linestyle ="dotted", color = "black")
    #converting into µs and adding all factors neglected in simulation
    plt.plot(hybrid_t , optimal_FID,color = hc.optimal_FID, label=r"Optimal hybrid FID with $c = " + f"{optimal_cut_off}" + r"$")
    plt.plot(hybrid_t , suboptimal_FID_1,color = hc.suboptimal_FID_1, label=r"Suboptimal hybrid FID with $c = " + f"{optimal_cut_off + suboptimal_factor1}" + r"$")
    plt.plot(hybrid_t , suboptimal_FID_2,color = hc.suboptimal_FID_2 ,label=r"Suboptimal hybrid FID with $c = " + f"{optimal_cut_off + suboptimal_factor2}" + r"$")
    #plt.plot(hybrid_t, QM_FID, color = "red", label=r"Gräßers data")
    #plt.plot(hybrid_t, classical_FID, color = "green", label=r"Classical data")
    #plt.plot(hybrid_t , optimal_FID,color = "blue")
    plt.plot(exp_time, exp_data, color = hc.exp_color, label=r"Experimental data")
    plt.xlabel( hc.time_label)
    plt.ylabel( hc.FID_label)
    plt.legend()
    plt.tight_layout()
    pdfname = "hybrid_" + pair_filename[:-4] + "_suboptimal" 
    plt.savefig("build/plots/hybrid_approach/" + pdfname + ".pdf")
    #check if plot is in master thesis folder and replace it if it is
    matching_paths = find_pdfs( pdfname)
    for path in matching_paths:
        print(f"\n\tReplacing plot in master thesis folder: {path}")
        plt.savefig(path )
    plt.close()


    return

def plot_cut_off_deviation(pair_filename):

    system_amount_name = extract_number(pair_filename[:10])
    particle_number_name = extract_number(pair_filename[10:])

    old_name = []
    if "_old" in pair_filename:
        old_name = "_old"
    else:
        old_name = ""

    #choose exp data
    if "[001]" in pair_filename:
        exp_data = exp_100
        exp_time = t_100
    elif "[011]" in pair_filename:
        exp_data = exp_110
        exp_time = t_110        
    elif "[111]" in pair_filename:
        exp_data = exp_111
        exp_time = t_111

    #calculate deviations for different cut offs
    deviation, std_dev = calc_different_cut_offs(exp_data, exp_time, pair_filename)

    #plotting
    size_label = hc.size_label

    plt.figure()
    plt.rc('axes', labelsize=size_label)
    #plt.figure(figsize=(10, 6))
    plt.errorbar(range(len(deviation)), deviation, yerr=std_dev, fmt='o', color = hc.hybrid_deviation, label='Deviation from experiment')# with std dev')
    plt.xlabel(hc.truncation_label)
    plt.xticks(np.arange(0,len(deviation),2,dtype=int))
    plt.ylabel(hc.deviation_3D_label_hybrid )
    plt.tight_layout()
    #plt.title(fr'Deviation for system with {system_amount_name} i.c. and {particle_number_name}$^3$ particles')
    plt.legend()
    
    #save figure
    pdfname = f"{pair_filename[:-4]}_cut_off_deviation"
    plt.savefig("build/plots/hybrid_approach/" + pdfname + ".pdf")
    #check if plot is in master thesis folder and replace it if it is
    matching_paths = find_pdfs( pdfname)
    for path in matching_paths:
        print(f"\n\tReplacing plot in master thesis folder: {path}")
        plt.savefig(path )
    plt.close()

    return deviation, std_dev


####include experimental data for FID
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










all_filenames = np.array([])

for file in os.listdir("build/daten"):
    all_filenames = np.append(all_filenames,file)



mask = [False] * len(all_filenames)
for i in range(0, len(all_filenames)):
    mask[i] = "txt" in all_filenames[i] and "pair" in all_filenames[i] and "interim" not in all_filenames[i] and "variance" not in all_filenames[i]
all_filenames = all_filenames[mask]

# just so its not empty and i can append
allready_checked_sys_amount = np.array([])
allready_checked_B_field = np.array([])

#this sorting is useless but i just copy and pasted it. maybe ist usefull some day
for file in all_filenames:
    
    system_amount = extract_number(file[:5])
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

    #make array to print a table to the terminal
    deviations = []
    std_devs = []

    for file in grouped_filenames:
        B_field_name =  "[" + file.split("[")[1].split("]")[0] + "]"



        
        plot_hybrid_approach(file)
        devs, standard_devs = plot_cut_off_deviation(file)
        if(devs is None or standard_devs is None):
            continue
        deviations.append(devs)
        std_devs.append(standard_devs)

    deviations = np.asarray(deviations)
    std_devs = np.asarray(std_devs)
