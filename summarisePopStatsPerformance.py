import sys
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Helper Function
def printNestedDict(dictionary):
    for effpopsize, models in dictionary.items():
        print(f"{effpopsize}:")
        for model_name, stats in models.items():
            print(f"  {model_name}:")
            for key, value in stats.items():
                print(f"   \t{key}:")
                for stat, values in value.items():
                    print(f"\t\t{stat}")
                    for pops, contents in values.items():
                        print(f"\t\t\t{pops}")
                        for key, items in contents.items():
                            print(f"\t\t\t\t{key}: {items[:4]}")

def grab_effpopsize_from_directory(directory):
    "Expects ./ and effpopsize is 2nd nested directory"
    effpopsize = directory.split('/')[2].split('_')[-1]
    return effpopsize

def intersect_dicts(dict1, dict2):
    # Placeholder values
    placeholder = {
        "wcFst": None,
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
        "wcFst": [dict1["wcFst"][i] for i in indices_to_retain],
        "pi": [dict1["pi"][i] if i < len(dict1["pi"]) else None for i in indices_to_retain], # protective indexing protects from unexpect lines
        "pi": [dict1["ehh"][i] if i < len(dict1["ehh"]) else None for i in indices_to_retain],
        "chromosome": [dict1["chromosome"][i] for i in indices_to_retain],
        "pos": [dict1["pos"][i] for i in indices_to_retain],
        "pi": [dict1["start_pos"][i] if i < len(dict1["start_pos"]) else None for i in indices_to_retain],
        "pi": [dict1["end_pos"][i] if i < len(dict1["end_pos"]) else None for i in indices_to_retain]
    }

    return intersected_dict

def captureGWASFilesIntoList(root_dir):
    # Input: base directory containing radinitio simulated directories
    # Loop through directories and append files to a list
    files_pattern = re.compile(r'(^nf-(vg|radseq)-(denovo|reference)_(pop(1|2|3)_pi|pop(1|2|3)_vs_pop(1|2|3)_wcFst)|^ri-master_(pops(1|2|3)_pi|pop(1|2|3)_vs_pop(1|2|3)_wcFst)).txt$')
    matched_files = []
    
    for dirpath, dirnames, filenames in os.walk(root_dir):
        dirnames[:] = [d for d in dirnames if d not in ['nf-radseq', 'work', 'vg']]
        for filename in filenames:
            if files_pattern.match(filename):
                full_path = os.path.join(dirpath, filename)
                matched_files.append(full_path)

    return matched_files

# Searches for specific files in directory
    # Outputs dictionary object with contents organized within effpopsize
    # Directory argument must not end with '/'
def grab_pop_statistics(file_path):
    stats = {'wcFst': [], 'pi': [], 'ehh': [], 'chromosome': [], 'pos': [], 'start_pos': [], 'end_pos': []}

    if os.path.isfile(file_path):
        with open(file_path, 'r') as file:
            for line in file:
                columns = line.strip().split('\t')
                
                # Check if columns length matches expected format
                if file_path.endswith('wcFst.txt') and len(columns) >= 3:
                    stats['chromosome'].append(columns[0])
                    stats['pos'].append(columns[1])
                    stats['wcFst'].append(float(columns[4]))

                elif file_path.endswith('pi.txt') and len(columns) >= 5:
                    stats['chromosome'].append(columns[0])
                    stats['pos'].append(columns[1])
                    stats['start_pos'].append(columns[1])
                    stats['end_pos'].append(columns[2])
                    stats['pi'].append(float(columns[3]))
                    stats['ehh'].append(float(columns[4]))

    # Sorting logic
    if file_path.endswith('wcFst.txt') and stats['wcFst']:
        zipped_lists = zip(stats['chromosome'], stats['pos'], stats['wcFst'])
        sorted_lists = sorted(zipped_lists, key=lambda x: (x[0], int(x[1])))
        stats['chromosome'], stats['pos'], stats['wcFst'] = zip(*sorted_lists)
    elif file_path.endswith('pi.txt') and stats['pi']:
        zipped_lists = zip(stats['chromosome'], stats['start_pos'], stats['end_pos'], stats['pi'], stats['ehh'])
        sorted_lists = sorted(zipped_lists, key=lambda x: (x[0], int(x[1]), int(x[2])))
        stats['chromosome'], stats['start_pos'], stats['end_pos'], stats['pi'], stats['ehh'] = zip(*sorted_lists)

    return stats

def extract_method(file):
    method = str(file).split('_')[0] # assuming located before first '_'
    return method

def extract_comparison(file):
    basename = file.split('/')[-1]
    if 'wcFst' in basename:
        comparison = basename.split('_')[1] + '_' + basename.split('_')[3]
    elif 'pi' in basename:
        comparison = basename.split('_')[1]
    else:
        raise ValueError(f'unrecognized command stat: {stat}')
    return comparison

def readGWASFiles(files):
    D = {}
    pi_seen = False  # Flag to check if 'pi' has been encountered

    for entry in files:
        effpopsize = grab_effpopsize_from_directory(entry)
        method = str(entry).split('/')[-1].split('.')[0].split('_')[0]
        category = "truth" if "ri-master" in entry else "model" # Determine if it's a truth or model data
        stat = str(entry).split('/')[-1].split('.')[0].split('_')[-1]  # Assuming 'wcFst', 'pi' or 'ehh' is the last part
        comparison = extract_comparison(entry)

        # Check if 'pi' has been seen before and if the current stat is 'pi' again
        #if pi_seen and stat == 'pi':
        #    stat = 'ehh'
        #elif stat == 'pi':
        #    pi_seen = True

        ref_data = grab_pop_statistics(entry) # dictionary with population statistics

        if 'wcFst' in entry:
            stat_data = grab_pop_statistics(entry.replace("_truth", ""))
            # Prune records based on overlap
            ref_data = intersect_dicts(ref_data, stat_data)

        # initialize dictionary if not already done so
        if effpopsize not in D:
            D[effpopsize] = {}
        if method not in D[effpopsize]:
            D[effpopsize][method] = {}
        if category not in D[effpopsize][method]:
            D[effpopsize][method][category] = {}
        if stat not in D[effpopsize][method][category]:
            D[effpopsize][method][category][stat] = {}
        if comparison not in D[effpopsize][method][category][stat]:
            D[effpopsize][method][category][stat][comparison] = {}

        D[effpopsize][method][category][stat][comparison] = ref_data
    return D

def extract_data(data_dict):
    plot_data = []
    metrics_of_interest = ['wcFst', 'pi', 'ehh']

    for effpopsize, highest_level in data_dict.items():
        for method, high_level in highest_level.items():
            for category, stat in high_level.items():
                for target_statistic in metrics_of_interest:  # Loop over each desired statistic
                    data = stat.get(target_statistic, {})  # Use .get() to get data or default to empty dict if not present
                    for comparison, pops in data.items():
                        for metric_stat, values in pops.items():
                            if metric_stat in metrics_of_interest:  # This check might be redundant now since we're looping over each metric
                                for value in values:
                                    plot_data.append((effpopsize, method, category, target_statistic, comparison, value))

    df = pd.DataFrame(plot_data, columns=["EffPopSize", "Method", "Category", "Statistic", "comparison", "Value"])
    df['Method'] = df['Method'].astype(str)
    df['EffPopSize'] = df['EffPopSize'].astype(int)
    df['Statistic'] = df['Statistic'].astype(str)
    df['combined_hue'] = df['Method'] + '_' + df['Category'].astype(str)

    return df
def flatten_dict(data):
    rows = []

    effpopsize_outer = list(data.keys())[0]
    first_level = data[effpopsize_outer]

    for effpopsize, second_level in first_level.items():
        for method, third_level in second_level.items():
            for category, fourth_level in third_level.items():
                for comparison, fifth_level in fourth_level.items():
                    for statistic, details in fifth_level.items():
                        for value in details:
                            num_entries = len(statistic.get('pos', []))
                            for i in range(num_entries):
                                row = {
                                    'EffPopSize': effpopsize,
                                    'Method': method,
                                    'Category': category,
                                    'Statistic': statistic,
                                    'Comparison': comparison,
                                    'wcFst': details.get('wcFst', [None] * num_entries)[i],
                                    'pi': details.get('pi', [None] * num_entries)[i],
                                    'ehh': details.get('ehh', [None] * num_entries)[i],
                                    'Chromosome': details.get('chromosome', [None] * num_entries)[i],
                                    'Pos': details.get('pos', [None] * num_entries)[i],
                                    'Start_Pos': details.get('start_pos', [None] * num_entries)[i],
                                    'End_Pos': details.get('end_pos', [None] * num_entries)[i]
                                }
                                rows.append(row)
    return rows
def flatten_dict_2(data):
    rows = []

    for effpopsize, first_level in data.items():
        for method, second_level in first_level.items():
            for category, third_level in second_level.items():
                for statistic, fourth_level in third_level.items():
                    for comparison, details in fourth_level.items():
                        num_entries = len(details.get('pos', []))  # Assume 'pos' is a good representative for the length
                        for i in range(num_entries):
                            try:
                                row = {
                                    'EffPopSize': effpopsize,
                                    'Method': method,
                                    'Category': category,
                                    'Statistic': statistic,
                                    'Comparison': comparison,
                                    'wcFst': details.get('wcFst', [None] * num_entries)[i],
                                    'pi': details.get('pi', [None] * num_entries)[i],
                                    'ehh': details.get('ehh', [None] * num_entries)[i],
                                    'Chromosome': details.get('chromosome', [None] * num_entries)[i],
                                    'Pos': details.get('pos', [None] * num_entries)[i],
                                    'Start_Pos': details.get('start_pos', [None] * num_entries)[i],
                                    'End_Pos': details.get('end_pos', [None] * num_entries)[i]
                                }
                                rows.append(row)
                            except IndexError:
                                # This will catch the error and continue with the next iteration
                                continue

    return rows

def drawBoxPlots(readGWASFiles_dictionary):
    metrics_of_interest = ['wcFst', 'pi', 'ehh']

    #print(df[df["Statistic"]].head())
    # Get unique values for all columns
    unique_values = {col: df[col].unique() for col in df.columns}

    # Print the unique values
    for col, values in unique_values.items():
        print(f"Unique values for {col}: {values}")
    
    for metric in metrics_of_interest:
        metric_df = df[df["Statistic"].str.contains(metric) & df["Value"].notna()]
        #print(df[df["Statistic"] == "ehh"].head())
        ehh_na_count = df[df["Statistic"] == "ehh"]["Value"].isna().sum()
        #print(f"Number of NaN values for ehh: {ehh_na_count}")
        if metric_df.empty:
            print(f"no data available for {metric}. Skipping plot.")
            continue
        
        unique_effpopsize_values = sorted(df['EffPopSize'].unique())
        if not unique_effpopsize_values:
            print(f"No unique 'EffPopSize' values found for {metric}. Skipping plot.")
            continue

        statistic = "fixation index" if 'wcFst' in metric else "nucleotide diversity (pi)" if 'pi' in metric else "extended haplotype homozygosity (ehh)"

        plt.figure(figsize=(12, 8))
        sns.boxplot(x="EffPopSize", y="Value", hue="combined_hue", data=metric_df, order=unique_effpopsize_values)
        
        plt.xlabel('effective population size')
        plt.ylabel(f'{statistic}')
        plt.title(f"Boxplot for {metric}")
        plt.legend(loc='upper right')
        plt.tight_layout()

        save_plotname = f'effpopsize_{metric}_boxplot.png'
        plt.savefig(save_plotname)
        plt.close()

def drawManhattanPlot(df):
    print(f'creating manhattan plots for {df}')
    
    # Loop through each unique combination of Statistic and Method
    for statistic in df['Statistic'].unique():
        for method in df['Method'].unique():
            
            subset = df[(df['Statistic'] == statistic) & (df['Method'] == method)]
            
            # Skip combinations without data
            if subset.empty:
                continue
            
            plt.figure(figsize=(12, 6))
            sns.scatterplot(data=subset, x=f'{position}', y='Value', hue='combined_hue', legend=False)
            plt.yscale('log')
            plt.ylabel('-log(Value)')
            plt.title(f'Manhattan Plot for {statistic} and {method}')
            
            # Save the plot
            filename = f"./{statistic}_{method}_manhattan_plot.png"
            plt.savefig(filename)
            print(f"Saved plot to {filename}")
            
            plt.close()  # Close the current figure to free up memory

if __name__ == "__main__":
    directory = sys.argv[1]

    files = captureGWASFilesIntoList(directory) # list of filenames
    #print(files)

    dict = readGWASFiles(files) # reads files and stores them in dictionary object

    #printNestedDict(dict) # view the contents of triple nested dict

    df = extract_data(dict)

    flattened_data = flatten_dict_2(dict)
    df = pd.DataFrame(flattened_data)
    df.head()

    #drawBoxPlots(df)

    #drawManhattanPlot(df)
    
