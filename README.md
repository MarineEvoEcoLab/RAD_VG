# Evaluating benefits of variant graph workflows in RADseq analysis

**Working Directory:** `/project/pi_jpuritz_uri_edu/radinitio`

## [radinitio.sh](./radinitio.sh) workflow script


### 1. Simulate Reads
Provided a reference sequence radinitio.sh will create a directory based on the file name and simulate RADseq reads using [radinitio](https://catchenlab.life.illinois.edu/radinitio/) software.

```bash
if [ ! -d "$OUTDIR" ] || [ -z "$(ls -A "$OUTDIR")" ]; then
    mkdir -p "$OUTDIR"
    echo "Directory created: $OUTDIR"
    # Simulating a ddRAD library:
    radinitio --simulate-all \
        --genome $GENOME \
        --chromosomes $CHROMOSOME_LIST \
        --out-dir $OUTDIR \
        --n-pops 4 --pop-eff-size 2500 --n-seq-indv 10 \
        --library-type ddRAD --enz PstI --enz2 MspI \
        --insert-mean 350 --insert-stdev 35 \
        --pcr-cycles 9 --coverage 20 --read-length 150

        bcftools norm -f $GENOME -Oz -o $TRUTH_NORM_VCF $TRUTH_VCF
        tabix -p vcf $TRUTH_NORM_VCF
else
    echo "Directory already exists: $OUTDIR"
fi
```

### 2. write dummy quality scores (.fa -> .fq)

convert .fa files outputed by radinitio to fq with [seqtk](https://github.com/lh3/seqtk)

```bash
for fa in "${RAD_READS_DIR}"/*.fa.gz; do
    fq=$(echo $fa | sed 's/.fa.gz/.fq.gz/')
    if [ -f "$fq" ] && [ -s "$fq" ]; then
        echo "$fq exists and is not empty."
    else
        #echo "$fq does not exist or is empty."
        module -q load seqtk/2023-05-11
        seqtk seq -F '@' $fa | gzip -c > $fq
    fi
done
```

### 3. write input samplesheet for nextflow workflow [radseq](https://github.com/Gabriel-A-Barrett/radseq)

```bash
paste -d',' <(for i in "${RAD_READS_DIR}"/*.1.fq.gz; do basename $i | cut -f1 -d'.' -; done)\
 <(ls "${RAD_READS_DIR}"/*.1.fq.gz) <(ls "${RAD_READS_DIR}"/*.2.fq.gz)\
 <(for i in "${RAD_READS_DIR}"/*.1.fq.gz; do echo "false"; done)\
 >> "${RAD_READS_DIR}/$DIRNAME.csv"
```

### 4. Run [radseq](https://github.com/Gabriel-A-Barrett/radseq) in reference and denovo mode

```bash
nextflow run Gabriel-A-Barrett/radseq -r dev --input "${RAD_READS_DIR}/$DIRNAME.csv" -c ~/config/unity.config --sequence_type 'PE'\
 --method 'reference' --genome $GENOME --outdir "${OUTDIR}/nf-radseq" -resume --save_reference_indices true
```

Denovo: Additional parameters **minreaddepth_withinindividual** and **minreaddepth_betweenindividual**
```bash
nextflow run Gabriel-A-Barrett/radseq -r dev --input "${RAD_READS_DIR}/$DIRNAME.csv" -c ~/config/unity.config --sequence_type 'PE'\
 --method 'denovo' --outdir "${OUTDIR}/nf-radseq" -resume --save_reference_indices 'true'
```

**Important** reference vcf from radinitio `${OUTDIR}/ref_loci_vars/ri_master.vcf.gz` is subsetted to regions found in merged alignment file from radseq.

```bash
bcftools query -R $NF_RADSEQ_REFERENCE_BED -f '%CHROM\t%POS\n' $TRUTH_NORM_VCF > $OUTDIR/test_accuracy/ri_master_norm_nf_ref_chrom_pos.txt
```


### 5. Pass radseq vcf to [vg_workflow](https://github.com/Gabriel-A-Barrett/nf-vg-pipeline) 

This workflow will construct a variant graph for alignment and genotyping.

```bash
nextflow run Gabriel-A-Barrett/nf-vg-pipeline --fasta ${GENOME} --fai "${REF_FAI}"\
 --vcf "${REF_VCF}" --tbi "${REF_TBI}"\
 --fq "${FQ}" -resume --output_mode 'symlink' -c ~/config/unity.config --outdir "${OUTDIR}/vg"
```
### 6. Calculate model sensitivity, precision, F1 score with [measureVcfAccuracy.py](measureVcfAccuracy.py)

Using pandas we can calculate true positives based on records found in both datasets, false positives based on records only found in only the model, and false negatives based on records only found in the reference using a pandas merge method.

```python
import sys
import pandas as pd

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

print(f"Sensitivity: {sensitivity:.3f}")
print(f"Precision: {precision:.3f}")
print(f"F1 Score: {f1_score:.3f}"
```

### TODO: How much parameter tweaking is enough 

1. apply multiple filters and use the reference with the best sensitivity and precision metrics
2. apply multiple parameters in denovo space and choose the best model performance. 
3. Calculate percent Fst error between reference and compare metrics obtained in model