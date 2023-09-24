import pandas as pd
import sys
import csv
import re
import numpy as np
import matplotlib.pyplot as plt
import os

#warnings.filterwarnings(action='ignore', category=RuntimeWarning)
"""
Second argument (true/false)
Switch 1 (outputBestPerformanceVCFName=true)
    find most accurate file within input directory
    prints to console vcf name with highest f1 score
Switch 2 (outputBestPerformanceVCFName=false)
    Grab most accurate files across effpopsize directories within test_accuracy
    Plot f1 scores for linear & vg across effpopsize
"""

def obtainReferenceFromInputDir(directory):
    parts = directory.split('/')
    if len(parts) > 1:
        reference = '/'.join(parts[:3]) + '/test_accuracy/ri_master_norm_dedup_rename_query.txt'
        return reference
    return None

def obtainMethodFromInputDir(directory):
    method = ''
    if ('nf-radseq-reference' in directory):
        method = 'nf-radseq-reference'
    elif ('nf-vg-reference' in directory):
        method = 'nf-vg-reference'
    elif ('nf-radseq-denovo' in directory):
        method = 'nf-radseq-denovo'
    elif ('nf-vg-denovo' in directory):
        method= 'nf-vg-denovo'
    return method

def calculate_accuracy(reference_input, model_input):
    # Load reference and model data
    reference_data = pd.read_csv(reference_input, sep='\t', names=['CHROM', 'POS', 'REF_ALLELE', 'ALT_ALLELE'])
    model_data = pd.read_csv(model_input, sep='\t', names=['CHROM', 'POS', 'REF_ALLELE', 'ALT_ALLELE'])

    merged = pd.merge(reference_data, model_data, on=['CHROM', 'POS'], how='outer', indicator=True)

    total = len(merged['CHROM'])
    true_positives = ((merged['_merge'] == 'both') & (merged['REF_ALLELE_x'] == merged['REF_ALLELE_y']) & (merged['ALT_ALLELE_x'] == merged['ALT_ALLELE_y'])).sum()
    false_positives = (merged['_merge'] == 'right_only').sum()
    false_negatives = (merged['_merge'] == 'left_only').sum()

    sensitivity = true_positives / (true_positives + false_negatives)
    precision = true_positives / (true_positives + false_positives)
    f1_score = 2 * (precision * sensitivity) / (precision + sensitivity)

    return f1_score, false_positives, precision, sensitivity, total

def find_highest_number_file(directory, reference, method):
    """
    searches for files w/ chrompos.txt
    calcs f1 scores
    appends to lits

    !!! Evaluates only files in the directory !!!
    """
    
    highest_number = float("-inf")
    highest_number_file = None
    highest_number_vcf = None
    data = []

    for filename in os.listdir(directory):
        if os.path.isfile(os.path.join(directory, filename)):
            if filename.endswith('chrompos.txt') and 'RADseqPopFilters' not in filename:
                effpopsize = grab_effpopsize_from_directory(directory)
                if 'reference' in directory:
                    f1, false_positives, precision, sensitivity, total = calculate_accuracy(reference, os.path.join(directory, filename))
                    
                    if 'nf-radseq/reference' in directory:
                        fmiss = filename.split('_')[-3]
                        mac = filename.split('_')[-2]
                        data.append({"filename": filename, 'effpopsize': effpopsize, "method": method, "f1": float('{:,.3f}'.format(f1)), "false_positives": false_positives, "fmiss": fmiss, "mac": mac, 'total': total})
                    else:
                        data.append({"filename": filename, 'effpopsize': effpopsize, "method": method, "f1": float('{:,.3f}'.format(f1)), "false_positives": false_positives, "fmiss": 'NA', "mac": 'NA','total': total})

                    if f1 > highest_number:
                        highest_number = f1
                        highest_number_file = filename
                        highest_number_vcf = highest_number_file.replace("_chrompos.txt",".vcf.gz")
                else:
                    if 'nf-radseq/denovo' in directory:
                        fmiss = float(filename.split('_')[-3])
                        mac = int(filename.split('_')[-2])
                        total = count_number_of_lines(os.path.join(directory, filename))
                        data.append({"filename": filename, 'effpopsize': effpopsize, "method": method, "f1": "NA", "false_positives": "NA", "fmiss": fmiss, "mac": mac, 'total': total})
                        #print(f'filename: {filename}\n\tfmiss: {fmiss}\n\tmac: {mac}')
                        if (fmiss == .95 and mac == 6):
                            highest_number_vcf = filename.replace("_chrompos.txt",".vcf.gz")


    return highest_number_vcf, data

def count_number_of_lines(file):
    with open(r"", 'r') as fp:
        line_count = len(fp.readlines())
    return line_count

def find_all_effpopsize_test_accuracy_dir_files(root_dir):
    """
    similar to find_highest_number_file but looks within nested directories
    """
    reference = None
    data = []  # Ensure the data list is initialized
    
    for dirpath, dirnames, filenames in os.walk(root_dir, topdown=False):
        effpopsize = grab_effpopsize_from_directory(dirpath)
        
        if 'test_accuracy' in dirpath:
            print("Directory Path: ", dirpath)
            print("Directories = ", dirnames)
            print("Files = ", filenames)
            print('-'*10)
            
            # First loop to find the reference
            for filename in filenames:
                if 'ri_master' in filename:
                    reference = os.path.join(dirpath, filename)
                    break  # Exit the loop once you've found the reference
            
            # Ensure a reference was found before processing other files
            if reference:
                # Second loop to process other files
                for filename in filenames:
                    if 'ri_master' not in filename:
                        prefix = grab_method_from_filename(filename)
                        if 'denovo' not in prefix:
                            f1, false_positives, precision, sensitivity, total = calculate_accuracy(reference, os.path.join(dirpath, filename))
                            data.append({
                                "filename": filename,
                                "effpopsize": effpopsize,
                                "prefix": prefix,
                                "total": total,
                                "f1": float('{:,.3f}'.format(f1)),
                                "precision": precision,
                                "sensitivity": sensitivity 
                            })
                        else:
                            data.append({
                                "filename": filename,
                                "effpopsize": effpopsize,
                                "prefix": prefix,
                                "total": total,
                                "f1": 'NA',
                                "precision": 'NA',
                                "sensitivity": 'NA' 
                            })
            reference = None # assuming resets value before re-assignment
    return data

def grab_method_from_filename(filename):
    prefix = filename.split('_')[0]
    return prefix

def grab_effpopsize_from_directory(directory):
    "Expects ./ and effpopsize is 2nd nested directory"
    effpopsize = directory.split('/')[2].split('_')[-1]
    return effpopsize

def variant_class_logger(directory, reference):
    """
    searches for variant class subset files
    """
    classes = ['_hom.chrompos', '_het.chrompos', '_snps.chrompos', '_sv.chrompos']
    predicted_filenames = ['nf-vg', 'nf-radseq']

    for class_suffix in classes:
        reference_filename = f'ri_master{class_suffix}'
        for filename in predicted_filenames:
            predicted_filename = f'{predicted_filename}{class_suffix}'

def writeCSV(data, directory):
    csv_file_name = 'vcfPerformance.csv'
    fields = ['filename','effpopsize', 'prefix', 'total_variant_amount', 'f1', 'precision', 'sensitivity']
    try: 
        # Write the data to the CSV file
        with open(csv_file_name, 'w', newline='') as csv_file:
            fieldnames = data[0].keys()
            writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

            writer.writeheader()
            for row in data:
                filename = row['filename']
                writer.writerow(row)

        #print(f'The data has been saved to {csv_file_name}')
    except IOError:
        print("I/O error")

def plot_f1_vs_fmiss_mac(data, directory):
    # Extract F1 scores, Fmiss, and Mac from the data
    f1_scores = [entry["f1"] for entry in data if "f1" in entry]
    fmiss_values = [entry["fmiss"] for entry in data if "fmiss" in entry]
    mac_values = [entry["mac"] for entry in data if "mac" in entry]
    
    # Sort data based on Fmiss and Mac values
    sorted_data = sorted(zip(fmiss_values, mac_values, f1_scores), key=lambda x: (x[0], x[1]))
    fmiss_values, mac_values, f1_scores = zip(*sorted_data)

    # Create a scatter plot
    plt.figure(figsize=(10, 6))
    plt.scatter(fmiss_values, mac_values, c=f1_scores, cmap='viridis', s=100, marker='o', edgecolors='k')
    plt.colorbar(label='F1 Score')
    
    # Set axis labels and title
    plt.xlabel('fmiss')
    plt.ylabel('mac')
    plt.title('F1 Score vs. Fraction missing genotypes and Minor allele count')

    # Sort the x-axis and y-axis
    plt.xticks(sorted(list(set(fmiss_values))))
    plt.yticks(sorted(list(set(mac_values))))

    # Show the plot
    plt.grid(True)
    
    save_filename = 'f1_vs_fmiss_mac.png'  # Replace with your desired filename

    # Save the plot as a PNG image
    plt.tight_layout()
    plt.savefig(os.path.join(directory, save_filename))

def plot_barchart_f1(data, directory):

    # Filter and sort the data for each prefix
    radseq_data = sorted([entry for entry in data if entry["prefix"] == "nf-radseq-reference"], key=lambda x: int(x["effpopsize"]))
    vg_data = sorted([entry for entry in data if entry["prefix"] == "nf-vg-reference"], key=lambda x: int(x["effpopsize"]))

    radseq_f1 = [entry["f1"] for entry in radseq_data]
    vg_f1 = [entry["f1"] for entry in vg_data]

    effpopsizes = [entry["effpopsize"] for entry in radseq_data]  # Assumes the effpopsize is same for both radseq and vg for the same index

    barWidth = 0.25
    fig, ax = plt.subplots(figsize=(12, 8))

    br1 = np.arange(len(radseq_f1))
    br2 = [x + barWidth for x in br1]

    ax.bar(br1, radseq_f1, color='b', width=barWidth, edgecolor='grey', label='radseq-reference')
    ax.bar(br2, vg_f1, color='r', width=barWidth, edgecolor='grey', label='vg-reference')

    ax.set_xlabel('Effective Population Size', fontweight='bold', fontsize=15)
    ax.set_ylabel('F1 Score', fontweight='bold', fontsize=15)
    ax.set_xticks([r + barWidth/2 for r in range(len(radseq_f1))])
    ax.set_xticklabels(effpopsizes)

    ax.legend()
    # Show the plot
    plt.grid(False)
    
    save_filename = 'f1_vs_method_vs_effpopsize.png'  # Replace with your desired filename

    # Save the plot as a PNG image
    plt.tight_layout()
    plt.savefig(os.path.join(directory, save_filename))


def main():

    if outputBestPerformanceVCFName:
        
        referenceFile = obtainReferenceFromInputDir(directory)
        method = obtainMethodFromInputDir(directory)
        best_vcf, data = find_highest_number_file(directory, referenceFile, method)

        print(best_vcf) # Return this value to console

        if method == 'nf-radseq-reference':      
            plot_f1_vs_fmiss_mac(data, directory)
       
    else:
        
        #referenceFile = obtainReferenceFromInputDir(directory)

        data = find_all_effpopsize_test_accuracy_dir_files(f'{directory}')
    
        # Summary DataFrame
        writeCSV(data, directory)
        
        # Visualize f1 scores for reference-linear & vg
        plot_barchart_f1(data, directory)

if __name__ == "__main__":
    
    # Get the input directory and output directory from command line arguments
    directory = sys.argv[1]
    outputBestPerformanceVCFName = True if len(sys.argv) > 2 and sys.argv[2].lower() == 'true' else False
    
    main()