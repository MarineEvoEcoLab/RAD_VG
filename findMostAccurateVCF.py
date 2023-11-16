import pandas as pd
import sys
import csv
import re
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns

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
    if directory.startswith('./'):
        reference = '/'.join(parts[:3]) + '/test_accuracy/ri-master_norm_dedup_rename_query.txt'
    else:
        reference = f'{directory}/' + '../../../../test_accuracy/ri-master_norm_dedup_rename_query.txt'
    return reference

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

    merged = pd.merge(reference_data, model_data, on=['CHROM', 'POS', 'REF_ALLELE', 'ALT_ALLELE'], how='outer', indicator=True)

    total = len(merged['CHROM'])
    true_positives = ((merged['_merge'] == 'both')).sum()
    false_positives = (merged['_merge'] == 'right_only').sum()
    false_negatives = (merged['_merge'] == 'left_only').sum()

    sensitivity = true_positives / (true_positives + false_negatives)
    precision = true_positives / (true_positives + false_positives)
    f1_score = 2 * (precision * sensitivity) / (precision + sensitivity)
    
    return f1_score, true_positives, false_positives, false_negatives, precision, sensitivity, total

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
    mask_files = ['RADseqPopFilters', 'dedup', 'impute']

    for filename in os.listdir(directory):
        if os.path.isfile(os.path.join(directory, filename)):
            if not filename.endswith('chrompos.txt') or 'RADseqPopFilters' in filename or 'dedup' in filename or 'impute' in filename or 'norm' in filename:
                continue
            else:
                file_path = os.path.join(directory, filename)
                if os.path.getsize(file_path) > 0:  # Check if the file is not empty
                    effpopsize, testing = grab_effpopsize_from_rep_directory(directory)
                    if 'reference' in directory:
                        f1, true_positives, false_positives, false_negatives, precision, sensitivity, total = calculate_accuracy(reference, file_path)
    
                        fmiss = filename.split('_')[-3]
                        mac = filename.split('_')[-2]
                        data.append({"filename": filename, 'effpopsize': effpopsize, "method": method, "f1": float('{:,.3f}'.format(f1)), "false_positives": false_positives, "fmiss": fmiss, "mac": mac, 'total': total})

                        if f1 > highest_number:
                            highest_number = f1
                            highest_number_file = filename
                            highest_number_vcf = highest_number_file.replace("_chrompos.txt",".vcf.gz")
                    elif 'denovo' in directory:
                        fmiss = float(filename.split('_')[-3])
                        mac = int(filename.split('_')[-2])
                        total = count_number_of_lines(file_path)
                        data.append({"filename": filename, 'effpopsize': effpopsize, "method": method, "f1": "NA", "false_positives": "NA", "fmiss": fmiss, "mac": mac, 'total': total})
                        #print(f'filename: {filename}\n\tfmiss: {fmiss}\n\tmac: {mac}')
                        if (fmiss == .35 and mac == 3):
                            highest_number_vcf = filename.replace("_chrompos.txt",".vcf.gz")
                    else: 
                        print

    return highest_number_vcf, data

def count_number_of_lines(file):
    with open(file, 'r') as fp:
        line_count = len(fp.readlines())
    return line_count

def find_all_effpopsize_test_accuracy_dir_files(root_dir):
    data = []  # Ensure the data list is initialized
    
    for dirpath, dirnames, filenames in os.walk(root_dir):
        # Check if the current directory is named 'test_accuracy'
        if os.path.basename(dirpath) == 'test_accuracy':
            print("Directory Path: ", dirpath)
            print("Directories = ", dirnames)
            print("Files = ", filenames)
            print('-'*10)
            
            effpopsize, repetition = grab_effpopsize_from_directory(dirpath)
            
            reference = None
            
            # First loop to find the reference
            for filename in filenames:
                if 'ri-master' in filename:
                    reference = os.path.join(dirpath, filename)
                    break  # Exit the loop once you've found the reference
            
            # Ensure a reference was found before processing other files
            if reference:
                # Second loop to process other files
                for filename in filenames:
                    if 'ri-master' not in filename:
                        prefix = grab_method_from_filename(filename)
                        if 'denovo' not in prefix:
                            f1, true_positives, false_positives, false_negatives, precision, sensitivity, total = calculate_accuracy(reference, os.path.join(dirpath, filename))
                            data.append({
                                "filename": filename,
                                "effpopsize": effpopsize,
                                "repetition": repetition,
                                "prefix": prefix,
                                "total": total,
                                'true_positives': true_positives,
                                'false_positives': false_positives,
                                'false_negatives': false_negatives,
                                "f1": float('{:,.3f}'.format(f1)),
                                "precision": precision,
                                "sensitivity": sensitivity 
                            })
                        else:
                            data.append({
                                "filename": filename,
                                "effpopsize": effpopsize,
                                "repetition": repetition,
                                "prefix": prefix,
                                "total": total,
                                'true_positives': true_positives,
                                'false_positives': false_positives,
                                'false_negatives': false_negatives,
                                "f1": 'NA',
                                "precision": 'NA',
                                "sensitivity": 'NA' 
                            })
            reference = None  # Reset the reference
        
    return data

def grab_method_from_filename(filename):
    prefix = filename.split('_')[0]
    return prefix

def grab_effpopsize_from_rep_directory(directory):
    # Get the current working directory
    current_directory = os.getcwd()
    # Get the parent directory
    parent_directory = os.path.dirname(current_directory)
    # Get the Lp_genome_Chr26_{type}
    grandparent_directory = os.path.dirname(parent_directory)

    return parent_directory, grandparent_directory

def grab_effpopsize_from_directory(directory):
    "Expects ./ and effpopsize is 2nd nested directory"
    effpopsize = directory.split('/')[-3]
    repetition = directory.split('/')[-2]
    return effpopsize, repetition

def variant_class_logger(directory, reference):
    """
    searches for variant class subset files
    """
    classes = ['_hom.chrompos', '_het.chrompos', '_snps.chrompos', '_sv.chrompos']
    predicted_filenames = ['nf-vg', 'nf-radseq']

    for class_suffix in classes:
        reference_filename = f'ri-master{class_suffix}'
        for filename in predicted_filenames:
            predicted_filename = f'{predicted_filename}{class_suffix}'

def writeCSV(data, directory):
    csv_file_name = directory + 'vcfPerformance.csv'
    fields = ['filename','effpopsize', 'prefix', 'total_variant_amount', 'false_positives', 'false_negatives', 'true_positives', 'f1', 'precision', 'sensitivity']
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

def plot_violinplot(data_list, directory):
    def extract_data_from_list(data_list):        
        plot_data = []
        for entry in data_list:
            effpopsize = entry['effpopsize']
            prefix = entry['prefix']
            filename = entry['filename']
            
            # We will be interested in both f1 and sensitivity
            metrics_of_interest = ['f1']
            
            for metric in metrics_of_interest:
                value = entry[metric]
                plot_data.append((effpopsize, prefix, filename, metric, value))
    
        return plot_data

    result = extract_data_from_list(data_list)
    df = pd.DataFrame(result, columns=['effpopsize', 'prefix', 'filename', 'metric', 'value'])

    # Convert 'effpopsize' to numeric
    df['effpopsize'] = pd.to_numeric(df['effpopsize'])

    # Get unique metrics and effective population sizes
    unique_metrics = df['metric'].unique()

    # Sort the DataFrame by 'effpopsize'
    df = df.sort_values(by='effpopsize')

    # Get unique effective population sizes after sorting
    unique_effpopsize_values = df['effpopsize'].unique()

    # Create a subplot for each metric
    fig, axs = plt.subplots(nrows=len(unique_metrics), figsize=(10, 6*len(unique_metrics)))

    for i, metric in enumerate(unique_metrics):
        ax = axs if len(unique_metrics) == 1 else axs[i]
        sns.violinplot(x='effpopsize', y='value', hue='prefix', 
                       data=df[df['metric'] == metric], order=unique_effpopsize_values, ax=ax)
        #ax.set_title(f'Violin plot of {metric} by Effective Population Size')
        ax.set_xlabel('Effective Population Size')
        ax.set_ylabel('F1 Score')
        save_filename = f'{metric}_violinplot.png'
        plt.tight_layout()
        plt.savefig(os.path.join(directory, save_filename))

def calculate_metrics(data):
    # Create a DataFrame from the provided data
    df = pd.DataFrame(data)

    # Convert 'effpopsize' to numeric (if not already)
    df['effpopsize'] = pd.to_numeric(df['effpopsize'])

    # Calculate metrics for all numeric columns
    numeric_columns = df.select_dtypes(include=['number']).columns
    metrics = ['mean', 'median', 'std']

    result = df.groupby(['filename', 'effpopsize'])[numeric_columns].agg(metrics).reset_index()
    
    # Rename the columns for clarity
    new_column_names = [f'{col}_{metric}' for col in numeric_columns for metric in metrics]
    result.columns = ['filename', 'effpopsize'] + new_column_names

    return result

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

        plot_violinplot(data, directory)

        print(calculate_metrics(data))

if __name__ == "__main__":
    
    # Get the input directory and output directory from command line arguments
    directory = sys.argv[1]
    
    outputBestPerformanceVCFName = True if len(sys.argv) > 2 and sys.argv[2].lower() == 'true' else False
    
    main()