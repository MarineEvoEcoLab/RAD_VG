import os
import sys
import gzip
import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re

def overlayed_barchart(data):
    df['effpopsize'] = pd.to_numeric(df['effpopsize'])

    positions = np.arange(len(df['effpopsize']))
    total_alignment = plt.bar(positions - width/2, df['total_alignments'], width, label='Total Alignments', zorder=1)
    properly_paired = plt.bar(positions + width/2, df['total_properly_paired'], width, label='Total Properly Paired', zorder=2)

    plt.xlabel('effpopsize')
    plt.ylabel('Number of Reads')
    plt.title('')
    plt.xticks(positions, df['effpopsize'])
    plt.legend()

    plt.tight_layout()


def violin_plot(df):
    print('violin plot')
    # Convert 'effpopsize' to numeric
    df['effpopsize'] = pd.to_numeric(df['effpopsize'])

    # Get unique metrics and effective population sizes
    unique_metrics = df['Type'].unique()

    # Sort the DataFrame by 'effpopsize'
    df = df.sort_values(by='effpopsize')

    # Get unique effective population sizes after sorting
    unique_effpopsize_values = df['effpopsize'].unique()

    # Define color palette
    color_palette = {'nf-radseq-reference': 'blue', 'nf-vg-reference': 'orange'}
    hue_order = ['nf-radseq-reference', 'nf-vg-reference']

    # Create a subplot for each metric
    fig, axs = plt.subplots(nrows=len(unique_metrics), figsize=(10, 6*len(unique_metrics)))

    for i, metric in enumerate(unique_metrics):
        ax = axs if len(unique_metrics) == 1 else axs[i]
        sns.boxplot(x='effpopsize', y='METRIC.F1_Score', hue='filename', hue_order=hue_order, 
                       data=df[df['Type'] == metric], order=unique_effpopsize_values, ax=ax, palette=color_palette)
        #ax.set_title(f'Violin plot of {metric} by Effective Population Size')
        ax.set_xlabel('Effective Population Size')
        ax.set_ylabel('F1 Score')
    save_filename = f'violinplot.png'
    plt.tight_layout()
    plt.savefig(os.path.join(directory, save_filename))
    plt.close()

def extract_genotype_scores(data, effpopsize, repetition):
    result_dict = {
        "filename": [],
        "effpopsize": [],
        "repetition": [],
    }

    # Extract the second string from the "value" field
    value_str = data.get("runInfo", [{}])[0].get("value", "")
    filename = value_str.split(" ")[-1].split('/')[-1]

    if "metrics" in data and len(data["metrics"]) > 0 and "data" in data["metrics"][0]:
        for item in data['metrics'][0]['data']:
            id_ = item['id']
            values = item['values']
            
            result_dict[id_] = values
        
        length = len(values)
        result_dict["filename"] = [filename] * length
        result_dict["effpopsize"] = [effpopsize] * length
        result_dict["repetition"] = [repetition] * length
        
        # Convert the dictionary to a pandas DataFrame
        df = pd.DataFrame(result_dict)
        # Filter the DataFrame based on the "Filter" column containing "PASS"
        df_filtered = df[df['Filter'] == 'PASS'].reset_index(drop=True)
        return df_filtered
    else:
        print("Metrics or data not found in data")
        return pd.DataFrame()

def calculate_summary_statistics(df):
    # Define the columns to calculate summary statistics for
    numeric_cols = df.select_dtypes(include='number').columns.tolist()
    
    # Group by 'filename', 'repetition', and 'effpopsize', and calculate summary statistics
    summary_stats = df.groupby(['filename', 'repetition', 'effpopsize', 'Type'])[numeric_cols].agg(['mean', 'std', 'min', 'max']).reset_index()

    return summary_stats

def read_json_file(file_path, human_readable=False):
    """
    Read a JSON file (either regular .json or .gz compressed) and return the parsed JSON data.
    Args:
        file_path (str): The path to the JSON file.
    Returns:
        dict: The parsed JSON data as a Python dictionary.
    """
    if file_path.endswith('.gz'):
        # Open the .gz file for reading in binary mode
        with gzip.open(file_path, 'rb') as gzipped_file:
            # Read the content of the file and decode it as UTF-8
            json_content = gzipped_file.read().decode('utf-8')
    else:
        # Open the file for reading
        with open(file_path, 'r') as json_file:
            # Read the content of the file
            json_content = json_file.read()

    # Parse the JSON content into Python objects
    json_data = json.loads(json_content)

    if human_readable:
        # If human_readable is True, format the JSON content with indentation
        json_content = json.dumps(json_data, indent=4)
    else:
        None

    return json_data, json_content

def grab_effpopsize_from_directory(directory_path):
    pattern = r"\b\d+\b"
    # Find all matches of the pattern in the directory_path
    matches = re.findall(pattern, directory_path)

    effpopsize = matches[0]
    repetition = matches[1]

    return effpopsize, repetition

def extract_alignment_stats(file_path, filename, effpopsize, repetition):
    results_dict = {
        'filename': filename,
        'effpopsize': effpopsize,
        'repetition': repetition
    }
    patterns = {
        "total_alignments": [r"(\d+) \+ \d+ properly paired", r"Total alignments: (\d+)"],
        "total_properly_paired": [r"(\d+) \+ \d+ mapped", r"Total properly paired: (\d+)"]
    }
    with open(file_path, 'r') as file:
        text = file.read()
    
    for key, regex_list in patterns.items():
        if not isinstance(regex_list, list):
            regex_list = [regex_list]
        
        value = None
        for regex in regex_list:
            match = re.search(regex, text)
            if match:
                value = int(match.group(1))
                break

        if value is not None:
            results_dict[key] = value
        else:
            results_dict[key] = None  # or you can choose to raise an exception

    return results_dict

def find_all_effpopsize_test_accuracy_dir_files(root_dir):
    """
    Parses through directory statistics on simulations and returns alignment, genotype, genome_scan data
    Args:
        root directory (str)
    Returns:
        multiple dict: 
    """
    
    data = []  # Ensure the data list is initialized
    filename_data_pairs = []
    test_accuracy = {}
    alignment_stats = {}
    genome_scan_stats = {}
    initial_depth = root_dir.rstrip(os.sep).count(os.sep)
    skip_dirs = ['vg','radseq',]

    for dirpath, dirnames, filenames in os.walk(root_dir):
        print(dirpath)
        depth = dirpath.rstrip(os.sep).count(os.sep) - initial_depth
        if depth == 2:
            effpopsize, repetition = grab_effpopsize_from_directory(dirpath)
            metric_data = {}  # Dictionary to store metric data
            metric_data['effpopsize'] = effpopsize
            metric_data['rep'] = repetition
        
        # Check for genotype performance directory
        if os.path.basename(dirpath) == 'test_accuracy':
            for filename in filenames:
                
                if '.metrics.json.gz' in filename:

                    json_data, json_content = read_json_file(os.path.join(dirpath, filename))
                    
                    df = extract_genotype_scores(json_data, effpopsize, repetition)
                    filename_data_pairs.append((filename, df))
                    print(f"effpopsize: {effpopsize}")
                    print(f"rep: {repetition}")
                    print(f"filename: {filename}")
                    print(df)
                    print('-'*10)
                    
                    metric_data['filename_data_pairs'] = filename_data_pairs
        
            data.append(metric_data)

        # Check for nf-core/radseq alignement stats directory
        elif os.path.basename(dirpath) == 'samtools_stat' or 'STAT' in os.path.basename(dirpath):
            for filename in filenames:
                if 'flagstat' in filename or '.stat' in filename:
                    metric_data['filename'] = filename
                    alignment_stats = extract_alignment_stats(os.path.join(dirpath, filename), filename, effpopsize, repetition)
    
        #elif os.path.basename(dirpath) == 'gwas':
            #for filename in filenames:


    # Combining genotype performance across replicates
    genotype_data_list = []

    for metric_data in data:
        if 'filename_data_pairs' in metric_data:
            filename_data_pairs = metric_data['filename_data_pairs']
            for filename, df in filename_data_pairs:
                genotype_data_list.append(df)

    if genotype_data_list:
        combined_df = pd.concat(genotype_data_list, ignore_index=True)
        print(combined_df)
        violin_plot(combined_df)

        calculate_metrics(combined_df,'./','low_gene_flow_genotype_performance.csv')
    else:
        print("No valid genotype data found.")


def calculate_metrics(data, directory, output_csv_file='genotypeSummary.csv', list_to_groupby = ['filename','effpopsize','Type']):
    # Create a DataFrame from the provided data
    df = pd.DataFrame(data)

    # Convert 'effpopsize' to numeric (if not already)
    df['effpopsize'] = pd.to_numeric(df['effpopsize'])

    # Calculate metrics for all numeric columns
    numeric_columns = df.select_dtypes(include=['number']).columns
    metrics = ['mean', 'median', 'std']

    result = df.groupby(list_to_groupby)[numeric_columns].agg(metrics).reset_index()
    
    # Rename the columns for clarity
    new_column_names = [f'{col}_{metric}' for col in numeric_columns for metric in metrics]
    result.columns = ['filename', 'effpopsize', 'Type'] + new_column_names
    # Save the result as a CSV file
    result.to_csv(os.path.join(directory, output_csv_file), index=False)

    return result

if __name__ == "__main__":

    directory = sys.argv[1]

    find_all_effpopsize_test_accuracy_dir_files(directory)