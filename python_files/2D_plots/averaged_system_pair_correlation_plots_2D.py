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
sys.path.append('/home/sander/Data/Masterarbeit/python_files/')

sys.path.append('/mnt/d/Dateien/Uni/Masterarbeit/python_files/')
import header_constants as hc



if os.path.exists("build") == False:
    os.mkdir("build")

if os.path.exists("build/plots2D") == False:
    os.mkdir("build/plots2D")

if os.path.exists("build/plots2D/averaged_pair_correlation") == False:
    os.mkdir("build/plots2D/averaged_pair_correlation")

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

    old_name = ""
    if("_old" in filename):
        old_name = "_old"

    #get pair correlations

    ##array for pair correlations
    #set to 0 for all correlations
    correlation_cut_off = 6
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

    if(correlation_cut_off > len(pair_correlations)):
        correlation_cut_off = len(pair_correlations)
    elif correlation_cut_off == 0:
        correlation_cut_off = len(pair_correlations)


    #calculate sum of pair correlations
    sum_of_all_correlations = pair_correlations[0]
    for i in range(1, len(pair_correlations)):
        sum_of_all_correlations = sum_of_all_correlations + pair_correlations[i]#np.add(sum_of_all_correlations, pair_correlations[i])

    print(f"\t\t{len(pair_correlations)} Different Pair Correlations found")
    
    C, err, t= np.genfromtxt("build/daten_2D/"+filename,skip_header= 1, unpack = True)


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



    plt.figure()
    plt.rc('axes', labelsize=size_label)
    plt.axhline(y = 0, xmin = 0, xmax = 1,linestyle ="dotted", color = "black")
    plt.plot(t , C, color ="black",label=r"Semi-classical FID")
    plt.plot(t[:len(sum_of_all_correlations)], sum_of_all_correlations,"--", color = "gold", label = r"$\sum m_c G_c^{xx}$")
    for i in range(0, correlation_cut_off):
        plt.plot(t[:len(pair_correlations[i])], pair_correlations[i], label=r"c = "  + str(i))
    plt.xlabel( hc.time_label_2D)
    plt.ylabel( hc.FID_label_2D)
    if(correlation_cut_off < 8 ):
        plt.legend()
    plt.legend()
    plt.tight_layout()
    pdfname = pair_filename[:-4]
    plt.savefig("build/plots2D/averaged_pair_correlation/" + pdfname + ".pdf")
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
            for pair_file in os.listdir("build/daten_2D"):
                if(pair_file.endswith("averaged_pair_Correlations.txt") == True and file[:34] in pair_file):
                    file_amount +=1        
                    print(f"\nFile {file} and {pair_file} opened\n")
                    plot_averaged_pair_correlation(file, pair_file)#, exp_100, t_100)
                    continue
        if(file.endswith("averaged_Correlation.txt") == True and "2D" in file and "[011]" in file and not "pair" in file):
            for pair_file in os.listdir("build/daten_2D"):
                if(pair_file.endswith("averaged_pair_Correlations.txt") == True and file[:34] in pair_file):
                    file_amount +=1        
                    print(f"\nFile {file} and {pair_file} opened\n")
                    plot_averaged_pair_correlation(file, pair_file)#, exp_110, t_110)
                    continue
        if(file.endswith("averaged_Correlation.txt") == True and "2D" in file and "[111]" in file and not "pair" in file):
            for pair_file in os.listdir("build/daten_2D"):
                if(pair_file.endswith("averaged_pair_Correlations.txt") == True and file[:34] in pair_file):
                    file_amount +=1        
                    print(f"\nFile {file} and {pair_file} opened\n")
                    plot_averaged_pair_correlation(file, pair_file)#, exp_111, t_111)
                    continue

        

        
print("Dateianzahl used for averaged pair correlation plots: ", file_amount)

