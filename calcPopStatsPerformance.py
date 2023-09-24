import sys
import os
import pandas as pd
import numpy as np
import glob
import re
import matplotlib.pyplot as plt
import seaborn as sns

# Searches for specific files in directory
    # Outputs dictionary object with contents organized within effpopsize
    # Directory argument must not end with '/'
def grab_pop_statistics(file_path):
    stats = {'pFst': [], 'pi': [], 'ehh': [], 'chromosome': [], 'pos': [], 'start_pos': [], 'end_pos': []}

    if os.path.isfile(file_path):
        with open(file_path, 'r') as file:
            for line in file:
                columns = line.strip().split('\t')
                stats['chromosome'].append(columns[0])
                stats['pos'].append(columns[1])
                
                if file_path.endswith('pFst.txt'):
                    stats['pFst'].append(float(columns[2]))
                elif file_path.endswith('pi.txt'):
                    stats['start_pos'].append(columns[1])
                    stats['end_pos'].append(columns[2])
                    stats['pi'].append(float(columns[3]))
                    stats['ehh'].append(float(columns[4]))
    # Sorting logic
    if file_path.endswith('pFst.txt'):
        zipped_lists = zip(stats['chromosome'], stats['pos'], stats['pFst'])
        sorted_lists = sorted(zipped_lists, key=lambda x: (x[0], int(x[1])))
        stats['chromosome'], stats['pos'], stats['pFst'] = zip(*sorted_lists)
    elif file_path.endswith('pi.txt'):
        zipped_lists = zip(stats['chromosome'], stats['start_pos'], stats['end_pos'], stats['pi'], stats['ehh'])
        sorted_lists = sorted(zipped_lists, key=lambda x: (x[0], int(x[1]), int(x[2])))
        stats['chromosome'], stats['start_pos'], stats['end_pos'], stats['pi'], stats['ehh'] = zip(*sorted_lists)

    return stats

def filter_stats_based_on_ref(ref_stats, target_stats):
    ref_positions = set(zip(ref_stats['chromosome'], ref_stats['pos']))
    
    # Get indices from target_stats that have positions in ref_positions
    indices_to_keep = [i for i, pos in enumerate(zip(target_stats['chromosome'], target_stats['pos'])) if pos in ref_positions]

    filtered_stats = {}
    for key in target_stats:
        filtered_stats[key] = [target_stats[key][i] for i in indices_to_keep if i < len(target_stats[key])]

    print("Length of ref_stats['chromosome']:", len(ref_stats['chromosome']))
    print("Length of target_stats['chromosome']:", len(target_stats['chromosome']))

    return filtered_stats

def flatten_dict(d):
    return {(chrom, pos): {"pFst": pFst, "pi": pi, "ehh": ehh, "start_pos": start, "end_pos": end} for chrom, pos, pFst, pi, ehh, start, end in zip(d["chromosome"], d["pos"], d["pFst"], d["pi"], d["ehh"], d["start_pos"], d["end_pos"])}

def intersect_dicts(dict1, dict2):
    # Placeholder values
    placeholder = {
        "pFst": None,
        "pi": None,
        "ehh": None,
        "start_pos": None,
        "end_pos": None
    }

    # Retain records from dict1 which also exist in dict2
    indices_to_retain = []
    for index, (chromosome, pos) in enumerate(zip(dict1['chromosome'], dict1['pos'])):
        if chromosome in dict2['chromosome'] and pos in dict2['pos']:
            indices_to_retain.append(index)
            
    # Populate result based on retained indices
    intersected_dict = {
        "pFst": [dict1["pFst"][i] for i in indices_to_retain],
        "pi": [dict1["pi"][i] if i < len(dict1["pi"]) else None for i in indices_to_retain], # protective indexing protects from unexpect lines
        "pi": [dict1["ehh"][i] if i < len(dict1["ehh"]) else None for i in indices_to_retain],
        "chromosome": [dict1["chromosome"][i] for i in indices_to_retain],
        "pos": [dict1["pos"][i] for i in indices_to_retain],
        "pi": [dict1["start_pos"][i] if i < len(dict1["start_pos"]) else None for i in indices_to_retain],
        "pi": [dict1["end_pos"][i] if i < len(dict1["end_pos"]) else None for i in indices_to_retain]
    }

    return intersected_dict

def grab_effpopsize_from_directory(directory):
    "Expects ./ and effpopsize is 2nd nested directory"
    effpopsize = directory.split('/')[2].split('_')[-1]
    return effpopsize

def collectReferenceSimulationsGWAS(root_dir):
    """
    looks within directory or *_truth*
    uses grab_pop_statistics
    """

    D_ref = {}
    D_stat = {}

    for directory_contents in glob.glob(root_dir + '/Lp_genome_*'):
        if os.path.isdir(directory_contents):
            effpopsize = grab_effpopsize_from_directory(directory_contents)
            
            ref_files = glob.glob(directory_contents + '/*_truth*')
            full_ref_files = glob.glob(directory_contents + '/ri_master*')
            all_ref_files = ref_files + full_ref_files
            for entry in all_ref_files:
                
                if effpopsize not in D_ref:
                    D_ref[effpopsize] = {}

                ref_data = grab_pop_statistics(entry)
                reference = str(entry).split('/')[-1].split('_')[0]

                if 'truth' in entry:
                    stat_filename = str(entry).replace("_truth", "")
                    
                    reference = str(entry).split('/')[-1].split('.')[0]
                    stat = stat_filename.split('/')[-1].split('.')[0]

                    if effpopsize not in D_stat:
                        D_stat[effpopsize] = {}

                    stat_data = grab_pop_statistics(stat_filename)
                    # If the filename contains pFst, prune records based on overlap
                    if 'pFst' in stat_filename:
                        new_stat_data = grab_pop_statistics(stat_filename)
                        stat_data = intersect_dicts(stat_data, ref_data)

                D_ref[effpopsize].setdefault(reference, {}).update(ref_data)
                D_stat[effpopsize].setdefault(stat, {}).update(stat_data)

    return D_ref, D_stat


# Sample D_stat structure
# D_stat = {
#     '1000': {
#         'pFst': {'chromosome': [...], 'pos': [...], ...},
#         'pi': {'chromosome': [...], 'pos': [...], ...},
#         ...
#     },
#     '2000': {...},
#     ...
# }


def splitDictionaryOnMethod(methods_dictionary):
    vg_dict = {}
    radseq_dict = {}
    for effpopsize, models in methods_dictionary.items():
        vg_inner_dict = {}
        radseq_inner_dict = {}
        for model_name, stats in models.items(): 
            if 'vg' in model_name:
                vg_inner_dict[model_name] = stats
            elif 'radseq' in model_name:
                radseq_inner_dict[model_name] = stats
        if vg_inner_dict:
            vg_dict[effpopsize] = vg_inner_dict
        if radseq_inner_dict:
            radseq_dict[effpopsize] = radseq_inner_dict
    return radseq_dict, vg_dict

def visualizePopStatsBoxPlot(ref_data, model_data):
    plot_data = []
    metrics = ['pFst', 'pi', 'ehh']

    def extract_data(data_dict):
        for effpopsize, methods in data_dict.items():
            for method, statistics in methods.items():
                for stat, values in statistics.items():
                    if stat in metrics:
                        for value in values:
                            plot_data.append((effpopsize, method, stat, value))

    # Extract data from both dictionaries and append to plot_data
    extract_data(ref_data)
    extract_data(model_data)

    df = pd.DataFrame(plot_data, columns=["EffPopSize", "Method", "Statistic", "Value"])
    # Convert the 'Method' and 'EffPopSize' columns to string and integer types, respectively
    df['Method'] = df['Method'].astype(str)
    df['EffPopSize'] = df['EffPopSize'].astype(int)
    
    for metric in metrics:
        # Filter the DataFrame for each metric
        metric_df = df[df["Statistic"] == metric]


        plt.figure(figsize=(12, 8))
        sns.boxplot(x="EffPopSize", y="Value", hue="Method", data=metric_df, order=sorted(df['EffPopSize'].unique()))
        
        # Customizing the plot
        plt.title(f"Boxplot for {metric}")
        plt.legend(loc='upper right')
        plt.tight_layout()

        # Save the plot
        save_plotname = f'effpopsize_{metric}_boxplot.png'
        plt.savefig(save_plotname)
        plt.close()  # Close the figure to free up memory

def visualizePopStatsDistribution(dict_ref, dict_radseq, dict_vg):
    print('building histograms')
    model= ['nf-radseq-reference', 'nf-vg-reference', 'nf-radseq-denovo', 'nf-vg-denovo']
    metrics = ['pFst', 'pi', 'ehh']
    colors = ['blue', 'green', 'red']
    for effpopsize, models in dict_radseq.items():
        for model_name, stats in models.items():
            radseq_model_name = model_name
            vg_model_name = model_name.replace('radseq','vg')
            reference_model_name = model_name.replace('nf-radseq-reference_','ri_master.')
            for key, values in stats.items():
                if key in metrics:
                    print(f"{effpopsize}:\n\t{model_name}:\n\t\t{key}: {values[:4]}")
                    radseq_values = dict_radseq[effpopsize][radseq_model_name][key]
                    vg_values = dict_vg[effpopsize][vg_model_name][key]
                    reference_values = dict_ref[effpopsize][reference_model_name][key]
                    for idx, metric in enumerate(metrics):
                        ax = axes[idx]
                            
                        if metric == key:
                            ax.hist(radseq_values, color=colors[0], alpha=0.6, label=f'{radseq_model_name} {metric}')
                            ax.hist(vg_values, color=colors[1], alpha=0.6, label=f'{vg_model_name} {metric}')
                            ax.hist(reference_values, color=colors[2], alpha=0.6, label=f'{reference_model_name} {metric}')
                                
                            ax.set_title(f'{metric} Distribution')
                            ax.set_xlabel(f'Value of {metric}')
                            ax.set_ylabel('Frequency')
                            ax.legend(loc='upper right')
                                
                            plt.tight_layout()
                            plt.subplots_adjust(top=0.85)

                            save_plotname = f'hist_{key}_{effpopsize}.png'

                            plt.savefig(save_plotname)

def plot_barchart_effpopsize_mse(data, directory):

    for stat in ['pFst_mse', 'pi_mse', 'ehh_mse']:
        
        # Extracting data for each statistic
        radseq_mse = [data[effpop]["nf-radseq-reference"][stat][0] for effpop in data]
        vg_mse = [data[effpop]["nf-vg-reference"][stat][0] for effpop in data]
        
        effpopsizes = list(data.keys())

        barWidth = 0.25
        fig, ax = plt.subplots(figsize=(12, 8))

        br1 = np.arange(len(radseq_mse))
        br2 = [x + barWidth for x in br1]

        ax.bar(br1, radseq_mse, color='b', width=barWidth, edgecolor='grey', label='radseq-reference')
        ax.bar(br2, vg_mse, color='r', width=barWidth, edgecolor='grey', label='vg-reference')

        ax.set_xlabel('Effective Population Size', fontweight='bold', fontsize=15)   
        ylabel = stat.split('_')[0]
        ax.set_ylabel(f'{ylabel} mean squared error (mse)', fontweight='bold', fontsize=15)
        ax.set_xticks([r + barWidth/2 for r in range(len(radseq_mse))])
        ax.set_xticklabels(effpopsizes)

        ax.legend()
        plt.grid(False)
        
        save_filename = f'{stat}_vs_method_vs_effpopsize.png' 

        plt.tight_layout()
        plt.savefig(os.path.join(directory, save_filename))
        plt.close()  # Close the figure to free up memory



def calculate_mse(actual_values, predicted_values):
    """
    Calculate the Mean Squared Error (MSE) between actual and predicted values.

    Parameters:
    actual_values (list or numpy.ndarray): List or array of actual values.
    predicted_values (list or numpy.ndarray): List or array of predicted values.

    Returns:
    float: Mean Squared Error (MSE) for the common elements between the two lists.
    """
    # Convert input to NumPy arrays for performance
    actual_values = np.array(actual_values)
    predicted_values = np.array(predicted_values)
    
        # Check if the input arrays have the same length
    if len(predicted_values) != len(actual_values):
        raise ValueError("Input arrays must have the same length.")
    # Find the common indices (elements that exist in both lists)
    #common_indices = np.intersect1d(np.arange(len(actual_values)), np.arange(len(predicted_values)))
    
    # Calculate MSE for the common elements
    squared_errors = (actual_values - predicted_values) ** 2
    mse = np.mean(squared_errors)
    
    return mse

def createModelPerformanceDictionary(actual_values, predicted_values):
    model_performance = {}
    for effpopsize, filename_dict in predicted_values.items():
        for mod_filename, stats in filename_dict.items():
            ref_filename, prefix = getLastElementFromSplit(mod_filename)
            #print(f'stat name used in actual_data: {ref_filename}')
            #print(prefix)
            if effpopsize not in model_performance:
                model_performance[effpopsize] = {}  # Initialize an empty dictionary for each effpopsize

            if prefix not in model_performance[effpopsize]:
                model_performance[effpopsize][prefix] = {'pFst_mse': [], 'pi_mse': [], 'ehh_mse': []}  # Initialize the 'fst_mse', 'pi_mse', 'ehh_mse' lists

            # Calculate and store the MSE for each statistic
            for key, value in stats.items():
                if key == 'pFst' or key == 'pi' or key == 'ehh':
                    
                    actual_data = actual_values.get(effpopsize, {}).get(ref_filename, {}).get(key, [])
                    predicted_data = predicted_values.get(effpopsize, {}).get(mod_filename, {}).get(key, [])

                    if actual_data and predicted_data and 'reference' in mod_filename:  # Check if both lists are non-empty
                        #print(f'effpopsize: {effpopsize}\tkey: {key}\n\treference: {ref_filename} (n={len(actual_data)})\tmodel: {mod_filename} (n={len(predicted_data)})')
                        mse = calculate_mse(actual_data, predicted_data)
                        model_performance[effpopsize][prefix][f'{key}_mse'].append(mse)
                        

    return model_performance

def getLastElementFromSplit(file):
    base_name, file_extension = os.path.splitext(file)
    base_name_list = base_name.split('_')  # Remove the leading '.' from the extension
    filename = base_name.split('.')[0]

    reference_name = re.sub(r"_(pi|pFst)$", r"_truth_\1", filename)

    prefix = base_name_list[0] 

    return reference_name, prefix    

def printDoubleNestedDict(dictionary):
    for effpopsize, models in dictionary.items():
        print(f"{effpopsize}:")
        for model_name, stats in models.items():
            print(f"  {model_name}:")
            for key, values in stats.items():
                print(f"    {key}: {values[:4]}")

def writePerformanceToCSV(data, output_filename):
    # Flatten the nested dictionary into a list of dictionaries
    append = True
    rows = []
    for effpopsize, filename in data.items():
        row = {'effpopsize': effpopsize}
        for model, stats in filename.items():
            for key, value in stats.items():
                # Check if the value is a list and convert it to a string without square brackets
                if isinstance(value, list):
                    row[f"{model}_{key}"] = ', '.join(map(str, value))
                else:
                    row[f"{model}_{key}"] = value
        rows.append(row)

    df = pd.DataFrame(rows)
    # Check if the file exists
    file_exists = os.path.isfile(output_filename)

    # If 'append' is True and the file exists, append to it; otherwise, create a new file
    if append and file_exists:
        # Append data to the existing file without writing headers
        df.to_csv(output_filename, mode='a', header=False, index=False)
    else:
        # Create a new file or overwrite the existing file
        df.to_csv(output_filename, index=False)


if __name__ == "__main__":
    
    print("Running calcPopStatsPerformance.py")
    directory = sys.argv[1]
    
    # Create a dictionary object for effpopsize -> workflow -> pop stats
    all_reference_stats, all_model_stats = collectReferenceSimulationsGWAS(directory)

    radseq_dict, vg_dict = splitDictionaryOnMethod(all_model_stats)

    printDoubleNestedDict(all_reference_stats)

    visualizePopStatsBoxPlot(all_reference_stats, all_model_stats)

    # Calculate mse
    performance = createModelPerformanceDictionary(all_reference_stats, all_model_stats)

    plot_barchart_effpopsize_mse(performance, directory)

    # write workflow mse scores
    writePerformanceToCSV(performance, 'finalSummary.csv')










