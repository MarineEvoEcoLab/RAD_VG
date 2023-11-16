# Evaluating performance between linear-reference and variant-graph workflows to measuring population summary statistics wcFst, pi, and ehh using RADseq.

**Working Directory:** `/project/pi_jpuritz_uri_edu/radinitio`

## [radinitio.sh](./radinitio.sh) workflow script

### Introduction

Despite advancements in generating genomic information it is still challenging to accurately classify genetic variation present, especially for complex regions. By considering population genetic variation at alignment and variant calling stages through graphical representation we can potentially reduce reference bias and improve detection of genetic variation. Using simulated restriction associated digest sequences from radinitio, we can better evaluate performance of variant graph to traditional linear-based methods. Provided a reference sequence, this software can create simulated RADseq, Run linear-based and variant-graph workflows, calculate genetic diversity and differientiation, and measure various workflow performance metrics. 

### Steps

- Write simulated RADseq based on reference sequences using radinitio
- Execute nf-core/radseq and Gabriel-A-Barrett/nf-vg-pipeline workflows with simulated individual RADseq in ${OUTDIR}/rad_reads based on reference and denovo methods. 
- Calculate workflow genotyping performance with hap.py
- Measure population summary statistics fixation index (pFst), nucleotide diversity (pi), and extended haplotype homozygosity (ehh) with vcflib
- Calculate summary statistics and visualize

### Prerequisites
- Linux/5.4.0-159-generic
- python/3.9.1
- radinitio/1.1.1
- nextflow/23.04.1.5866
- apptainer/1.1.5+py3.8.12

### Installation/Setup

To install software simply clone this repository onto your machine:
```
git clone <this repository>
```

Provide executable permissions to scripts:

```
chmod +x *.sh & chmod +x *.py
```

### Data
- Example data is provided within `./atlantic_hsc/Lp_genome_Chr26.fasta`


## Method Overview 
#### **Parameters**
```bash
command: bash radinitio.sh -g ./atlantic_hsc/Lp_genome_Chr26.fasta -n 2 -i 10 -o out 
    -g | --genome : reference sequence
    -n | --number_of_populations : number of simulated populations
    -i | --number_of_individuals : number of simulated individuals per-population 
    -o | --outdir : output directory
    -h | --help : display help message
     
```
Quick overview: 
1. Simulate reads
2. Create a bed file with regions from file: `./ref_loci_vars/reference_rad_loci.fa.gz` headers
3. Convert fa to fq files
4. Write samplesheet for linear pipeline
5. Run linear pipeline using `reference` method
5. Run graph pipeline from linear output
6. Scan the genome and calculate Fst, pi, ehh



### 1. Simulate Reads

#### 1.a modify radinitio to include and run simple_msp_stepping_stone_model in `__init__.py`

```python
#
# Make msprime stepping stone models
def simple_msp_stepping_stone_model(n_seq_indv, pop_eff_size, n_pops, mutation_rate=7e-8, pop_immigration_rate=0.001, growth_rate=0.0, rate=0.8):
    # Msprime parameters
    # Compute population sizes using list comprehension
    population_sizes = [pop_eff_size * (rate ** i) for i in range(n_pops)]

    pops = [
        msprime.PopulationConfiguration(
            sample_size = 2 * n_seq_indv,
            initial_size = population_sizes[i],
            growth_rate = growth_rate,)
        for i in range(n_pops) 
    ]
    
    msprime_simulate_args = {
        'mutation_rate' : mutation_rate,
        'population_configurations' : pops }
    
    # Create the migration matrix for stepping stone model.
    if n_pops > 1:
        migration_matrix = [ [0] * n_pops for _ in range(n_pops) ]
        for j in range(n_pops):
            if j > 0:
                migration_matrix[j][j-1] = pop_immigration_rate
            if j < n_pops - 1:
                migration_matrix[j][j+1] = pop_immigration_rate
        msprime_simulate_args['migration_matrix'] = migration_matrix # append to simulate args
    
    return msprime_simulate_args
```

and edited `__main__.py` to raise stepping stone model 

```python
# Define msprime population options
msprime_simulate_args = ri.simple_msp_stepping_stone_model(args.n_seq_indv, args.pop_eff_size, args.n_pops)
```

#### 1.b Run radinitio

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
        --n-pops 3 --pop-eff-size 2500 --n-seq-indv 30 \
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

### 3. subset truth regions 

```bash
singularity exec ${SIF_DIR}/htslib%3A1.18--h81da01d_0 sh -c "bgzip -d -c $TRUTH_FASTA | grep '>' | cut -f 2 -d'=' | sed 's/[:-]/\t/g' > ${TRUTH_REGIONS}" # create bed file
singularity exec ${SIF_DIR}/bcftools%3A1.17--haef29d1_0 bcftools view -Oz -o $TRUTH_BGZIP_VCF $TRUTH_VCF
singularity exec ${SIF_DIR}/htslib%3A1.18--h81da01d_0 tabix -p vcf -f $TRUTH_BGZIP_VCF # write .tbi
singularity exec ${SIF_DIR}/bcftools%3A1.17--haef29d1_0 bcftools view -R ${TRUTH_REGIONS} -Oz -o $TRUTH_SUB_VCF $TRUTH_BGZIP_VCF
singularity exec ${SIF_DIR}/htslib%3A1.18--h81da01d_0 tabix -p vcf -f $TRUTH_SUB_VCF # write .tbi
singularity exec ${SIF_DIR}/bcftools%3A1.17--haef29d1_0 bcftools norm -a -m -any -f $GENOME $TRUTH_SUB_VCF > $TRUTH_SUB_NORM_VCF
singularity exec ${SIF_DIR}/bcftools%3A1.17--haef29d1_0 bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $TRUTH_SUB_NORM_VCF > "${OUTDIR}/test_accuracy/ri-master_norm_dedup_rename_query.txt"
```


### 3. write input samplesheet for nextflow workflow [radseq](https://github.com/Gabriel-A-Barrett/radseq)

```bash
paste -d',' <(for i in "${RAD_READS_DIR}"/*.1.fq.gz; do basename $i | cut -f1 -d'.' -; done)\
 <(ls "${RAD_READS_DIR}"/*.1.fq.gz) <(ls "${RAD_READS_DIR}"/*.2.fq.gz)\
 <(for i in "${RAD_READS_DIR}"/*.1.fq.gz; do echo "false"; done)\
 >> "${RAD_READS_DIR}/$DIRNAME.csv"
```

### 4. Run [radseq](https://github.com/Gabriel-A-Barrett/radseq) in reference and denovo mode

expect pipelines are in your nextflow home bin folder

**Reference**
```bash
nextflow run Gabriel-A-Barrett/radseq -r dev --input "${RAD_READS_DIR}/$DIRNAME.csv" -c ~/config/unity.config --sequence_type 'PE'\
 --method 'reference' --genome $GENOME --outdir "${OUTDIR}/nf-radseq" -resume --save_reference_indices true
```

**Denovo**: Additional parameters **minreaddepth_withinindividual** and **minreaddepth_betweenindividual**

```bash
nextflow run Gabriel-A-Barrett/radseq -r dev --input "${RAD_READS_DIR}/$DIRNAME.csv" -c ~/config/unity.config --sequence_type 'PE'\
 --method 'denovo' --outdir "${OUTDIR}/nf-radseq" -resume --save_reference_indices 'true'
```

### 5. Build variant graph [vg_workflow](https://github.com/Gabriel-A-Barrett/nf-vg-pipeline) 

Uses [freebayes](https://github.com/freebayes/freebayes) and [vg pack](https://github.com/vgteam/vg) to call genotypes

```bash
nextflow run Gabriel-A-Barrett/nf-vg-pipeline --fasta ${GENOME} --fai "${REF_FAI}"\
 --vcf "${REF_VCF}" --tbi "${REF_TBI}"\
 --fq "${FQ}" -resume --output_mode 'symlink' -c ~/config/unity.config --outdir "${OUTDIR}/vg"
```
### 6. Evaluate model performance with [hap.py](https://github.com/Illumina/hap.py)
```bash
singularity exec ${SIF_DIR}/hap.py%3A0.3.15--py27hcb73b3d_0 hap.py <truth vcf> <model vcf> -r <reference fasta> -o <outdir>/test_accuracy/nf-<pipeline>-reference
```

### 7. Measure genome-wide Fst, pi, and ehh with [vcflib]()

```bash
extract_population_indices_and_calculate() {
    # Check if at least 3 arguments are provided
    if [ "$#" -lt 3 ]; then
        echo "Usage: extract_population_indices_and_calculate <path_to_vcf> <output_directory> <prefix> [ <number_of_populations> <number_of_individuals_per_population> ]"
        exit 1
    fi

    VCF="$1"
    OUTDIR="$2"
    prefix="$3"
    NUM_POPULATIONS="${4:-4}"
    INDV_PER_POP="${5:-30}"

    # Ensure the VCF file exists
    if [ ! -f "${VCF}" ]; then
        echo "VCF file not found!"
        exit 1
    fi

    # Check if OUTDIR exists, if not create it
    if [ ! -d "${OUTDIR}/gwas" ]; then
        echo "Output directory does not exist. Creating it..."
        mkdir -p "${OUTDIR}/gwas"
    fi

    # Get list of individuals with line numbers
    singularity exec ${SIF_DIR}/bcftools%3A1.17--haef29d1_0 bcftools query -l "${VCF}" | awk -v OFS='\t' '{print NR, $0}' | sort -k2 > "${VCF}_INDVS"

    # Extract individual names for each population and store in an array
    declare -a pop_indices_array
    for i in $(seq 1 $NUM_POPULATIONS); do
        start=$((($i-1) * $INDV_PER_POP + 1))
        end=$(($i * $INDV_PER_POP))
        pop_names=$(awk -v start="$start" -v end="$end" 'NR < start || NR > end { next } { printf("%s%s", sep, $1); sep="," } END { printf("\n") }' "${VCF}_INDVS")
        pop_indices_array[$i]="$pop_names"
    done

    # Run the command for every unique combination of populations, excluding self-comparisons
    for i in $(seq 1 $NUM_POPULATIONS); do
        for j in $(seq $i $NUM_POPULATIONS); do
            if [ $i -eq $j ]; then
                continue
            fi
            pops="pop${i}_vs_pop${j}"
            echo "Processing combination: $pops"

            if [ ! -f "$OUTDIR/${prefix}_${pops}_wcFst.txt" ]; then
                singularity exec ${SIF_DIR}/vcflib%3A1.0.3--hecb563c_1 wcFst --target "${pop_indices_array[$i]}" --background "${pop_indices_array[$j]}" -y GT --file "${VCF}" > "$OUTDIR/gwas/${prefix}_${pops}_wcFst.txt"
            else
                echo "$OUTDIR/${prefix}_${pops}_wcFst.txt already exists. Skipping."
            fi
        done
    done

```

### Output Directory Structure

- **test_accuracy** folder contains performance from *freebayes linear*, *freebayes graph alignment*, and *vg call*
- Genome scan files from 3 variant callers get deposited into **gwas** folder


```
├── {outdir}
│   ├── {effpopsize}
│   │   ├── {repetition}
│   │   │   ├── gwas
│   │   │   ├── msprime_vcfs
│   │   │   ├── nf-radseq
│   │   │   │   ├── fastp
│   │   │   │   ├── fastqc
│   │   │   │   ├── pipeline_info
│   │   │   │   └── reference
│   │   │   │       ├── bwa-mem2
│   │   │   │       │   ├── intervals
│   │   │   │       │   │   ├── bedops_merge
│   │   │   │       │   │   ├── bedtools_bamtobed
│   │   │   │       │   │   └── bedtools_makewindows
│   │   │   │       │   ├── samtools_index
│   │   │   │       │   ├── samtools_merge
│   │   │   │       │   └── samtools_stats
│   │   │   │       ├── multiqc
│   │   │   │       │   └── multiqc_data
│   │   │   │       ├── reference
│   │   │   │       │   └── bwa-mem2
│   │   │   │       │       └── index
│   │   │   │       ├── samtools
│   │   │   │       │   └── index
│   │   │   │       └── variant_calling
│   │   │   │           ├── filter
│   │   │   │           └── intervals
│   │   │   ├── rad_alleles
│   │   │   ├── rad_reads
│   │   │   ├── ref_loci_vars
│   │   │   ├── test_accuracy
│   │   │   └── vg
│   │   │       └── reference
│   │   │           ├── BCFTOOLS
│   │   │           │   └── MERGE
│   │   │           ├── STATS
│   │   │           ├── TABIX
│   │   │           │   └── TABIX
│   │   │           └── VG
│   │   │               └── PATH
|   |   |
```
