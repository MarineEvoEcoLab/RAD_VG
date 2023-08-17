import sys
import pandas as pd
from sklearn.metrics import roc_curve


reference_input=sys.argv[1]
model_input=sys.argv[2]

# Load reference and model positions
reference_positions = pd.read_csv(reference_input, sep='\t', names=['CHROM', 'POS'])
model_positions = pd.read_csv(model_input, sep='\t', names=['CHROM', 'POS'])

# Merge dataframes on 'CHROM' and 'POS'
merged = pd.merge(reference_positions, model_positions, on=['CHROM', 'POS'], how='outer', indicator=True)

# Calculate accuracy metrics
true_positives = (merged['_merge'] == 'both').sum()
false_positives = (merged['_merge'] == 'right_only').sum()
false_negatives = (merged['_merge'] == 'left_only').sum()

print(f"true_positives: {true_positives:.2f}")
print(f"false_negatives: {false_positives:.2f}")

sensitivity = true_positives / (true_positives + false_negatives)
precision = true_positives / (true_positives + false_positives)
f1_score = 2 * (precision * sensitivity) / (precision + sensitivity)


reference_positions['SCORE'] = reference_positions['CHROM'] + '_' + reference_positions['POS']
model_positions['SCORE'] = model_positions['CHROM'] + '_' + model_positions['POS']

fpr, tpr, thresholds = roc_curve(reference_positions['SCORE'],model_positions['SCORE'])

print(f"Sensitivity: {sensitivity:.3f}")
print(f"Precision: {precision:.3f}")
print(f"F1 Score: {f1_score:.3f}"