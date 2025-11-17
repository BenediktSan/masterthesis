import numpy as np
import matplotlib.pyplot as plt
#import uncertainties as unc
#import uncertainties.unumpy as unp
#from uncertainties import ufloat
from scipy.optimize import curve_fit
import scipy.constants as const
#import sympy
import os
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


def plot_averaged_pair_correlation(filename, pair_filename):#,FID_exp, t_exp):

    #get pair correlations

    ##array for pair correlations
    #set to 0 for all correlations
    correlation_cut_off = 6
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

    if(correlation_cut_off > len(pair_correlations)):
        correlation_cut_off = len(pair_correlations)
    elif correlation_cut_off == 0:
        correlation_cut_off = len(pair_correlations)


    #calculate sum of pair correlations
    sum_of_all_correlations = pair_correlations[0]
    for i in range(1, len(pair_correlations)):
        sum_of_all_correlations = sum_of_all_correlations + pair_correlations[i]#np.add(sum_of_all_correlations, pair_correlations[i])

    print(f"\t\t{len(pair_correlations)} Different Pair Correlations found")
    
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

    #normalizing
    norm_C = C[0]
    #C = C/norm_C

    norm = sum_of_all_correlations[0]
    #sum_of_all_correlations = sum_of_all_correlations / norm
    #pair_correlations = pair_correlations / norm

    print("\t\tNormalization factor for pair correlations: ", norm)
    print("\t\tNormalization factor for C: ", norm_C)
    #print(pair_correlations[0][0])
    #print(sum_of_all_correlations[0])


    size_label = hc.size_label

    #print(f"\t\tsize of t_tilde {t_tilde.size} and size of pair_corrs {pair_correlations[0].size} {pair_correlations[1].size}")


    plt.figure()
    plt.rc('axes', labelsize=size_label)
    plt.axhline(y = 0, xmin = 0, xmax = 1,linestyle ="dotted", color = "black")
    #converting into µs and adding all factors neglected in simulation
    #plt.plot(t_tilde , C,"--", color = "gold",label=r"$C_x$")
    #plt.fill_between(t_tilde, C + err, C - err, color = "green", alpha = 0.3, label ="estimated MC error")
    #plt.plot(t_tilde_2 , C, label=r"C3", alpha = 0.5)
    #plt.plot(t_exp, FID_exp, color ="purple", label=r"experimental data")
    plt.plot([] , [], color =hc.FID_color_pair_corr,label=r"Semi-classical FID")
    plt.plot([] , [], color = hc.pair_corr_sum_color, label = r"$\sum m_c G_c^{xx}(t)$")
    for i in range(0, correlation_cut_off):
        plt.plot(t_tilde[:len(pair_correlations[i])], pair_correlations[i], label=r"$c$="+str(i+1))
    plt.plot(t_tilde , C, color =hc.FID_color_pair_corr)
    plt.plot(t_tilde[:len(sum_of_all_correlations)], sum_of_all_correlations,"--", color = hc.pair_corr_sum_color)
    plt.xlabel( hc.time_label)
    plt.ylabel( hc.FID_label )
    if(correlation_cut_off < 8 ):
        plt.legend()
    plt.legend()
    plt.tight_layout()
    pdfname = pair_filename[:-4]
    plt.savefig("build/plots/averaged_pair_correlation/" + pdfname + ".pdf")
    #check if plot is in master thesis folder and replace it if it is
    matching_paths = find_pdfs( pdfname)
    for path in matching_paths:
        print(f"\n\tReplacing plot in master thesis folder: {path}")
        plt.savefig(path )
    plt.close()

    return



'''
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


'''


file_amount = 0

for file in os.listdir("build/daten"):
    if file.endswith(".txt"):
        #print(file)
        if(file.endswith("averaged_Correlation.txt") == True and "[001]" in file and not "pair" in file):
            for pair_file in os.listdir("build/daten"):
                if(pair_file.endswith("averaged_pair_Correlations.txt") == True and file[:42] in pair_file):
                    file_amount +=1        
                    print(222)
                    print(f"\nFile {file} and {pair_file} opened\n")
                    plot_averaged_pair_correlation(file, pair_file)#, exp_100, t_100)
                    continue
        if(file.endswith("averaged_Correlation.txt") == True and "[011]" in file and not "pair" in file):
            for pair_file in os.listdir("build/daten"):
                if(pair_file.endswith("averaged_pair_Correlations.txt") == True and file[:42] in pair_file):
                    file_amount +=1        
                    print(f"\nFile {file} and {pair_file} opened\n")
                    plot_averaged_pair_correlation(file, pair_file)#, exp_110, t_110)
                    continue
        if(file.endswith("averaged_Correlation.txt") == True and "[111]" in file and not "pair" in file):
            for pair_file in os.listdir("build/daten"):
                if(pair_file.endswith("averaged_pair_Correlations.txt") == True and file[:42] in pair_file):
                    file_amount +=1        
                    print(f"\nFile {file} and {pair_file} opened\n")
                    plot_averaged_pair_correlation(file, pair_file)#, exp_111, t_111)
                    continue

        

        
print("Dateianzahl used for averaged pair correlation plots: ", file_amount)

