# Comparing genotype and estimating population statistics performance between linear-reference and variant-graph workflows in RADseq.

**Working Directory:** `/project/pi_jpuritz_uri_edu/radinitio`

## [radinitio.sh](./radinitio.sh) workflow script

### Introduction

Despite advancements in generating genomic information it is still challenging to accurately classify genetic variation present, especially for complex regions. By considering population genetic variation at alignment and variant calling stages through graphical representation we can potentially reduce reference bias and improve detection of genetic variation. Using simulated restriction associated digest sequences from radinitio, we can better evaluate performance of variant graph to traditional linear-based methods. Provided a reference sequence, this software can create simulated RADseq, Run linear-based and variant-graph workflows, calculate genetic diversity and differientiation, and measure various workflow performance metrics. 

### Objectives

- Write simulated RADseq based on reference sequences using radinitio
- Execute nf-core/radseq and Gabriel-A-Barrett/nf-vg-pipeline workflows with simulated individual RADseq in ${OUTDIR}/rad_reads based on reference and denovo methods. 
- Calculate workflow genotyping performance: f1 score, precision, ect. with in-house python scripts
- Measure population summary statistics fixation index (pFst), nucleotide diversity (pi), and extended haplotype homozygosity (ehh) using vcflib
- Calculate summary stastic performance mean squared error and distribution of data. 

### Prerequisites
- Linux/5.4.0-159-generic
- python/3.9.1
- radinitio/1.1.1
- nextflow/23.04.1.5866
- apptainer/1.1.5+py3.8.12

### Installation/Setup

### Data
- Provide a ready to use reference sequences for testing? 


## Method Overview

#### Parameters
```
radinitio.sh
    argument $1: reference sequence
    argument $2: effective population size for simulations
```

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


### Results and Outputs

#### Examples

### Troubleshooting

### Contribution & Feedback

### Citations & Acknowledgments

### License

### Updates & Version History

### Contact Information


### TODO: How much parameter tweaking is enough 

1. Create a bash function that checks to see if the singularity images are downloaded
1. output markdown dataframe describing genotype performance and population summary statistics
2. 