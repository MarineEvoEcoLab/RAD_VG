#!/usr/bin/env bash
#SBATCH -c 20
#SBATCH --mem=10G
#SBATCH -p uri-cpu
#SBATCH --time=72:00:00
#SBATCH --out=./err_out/%x_%A_%a.out
#SBATCH --err=./err_out/%x_%A_%a.err
set -e # Exit on error

usage() {
    cat <<EOF
Usage: ${0##*/} [OPTIONS]

Options:
  -g, --genome                      Reference fasta sequence to simulate off of.
  -n, --number_of_populations       Number of populations to simulate.
  -i, --number_of_individuals       Number of individuals within population(s).
  -o, --outdir                      Output directory; suffix gets appended to reference genome name provided.
  -h, --help                        Display this help and exit.

Tasks are:
  1. Simulate rad sequences
     1.1 Calculate summary statistics
  2. Run radseq
     2.1 Calculate summary statistics
  3. Run vg
     3.1 Calculate summary statistics
EOF
}

# Check if the user asked for the usage information with -h or --help
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    usage
    exit 0
fi

# Check if the minimum number of arguments is provided
if [ "$#" -lt 0 ]; then
    usage
    exit 1
fi

module -q load apptainer/1.1.5+py3.8.12

eval SIF_DIR="/project/pi_jpuritz_uri_edu/radinitio/sif"
export NXF_SINGULARITY_CACHEDIR="/home/gabriel_barrett_uri_edu/nxf-singularity-cache-dir"

nextflow pull Gabriel-A-Barrett/radseq
nextflow pull Gabriel-A-Barrett/nf-vg-pipeline

vcf_rmdup_fixgt_phase() {
    # capture final output: phased vcf w/ tail
    # 
    local VCF=$1
    local RM_DUP_VCF=$(echo $VCF | sed 's/.vcf.gz/_dedup.vcf/') # append _dedup
    local BIALLELIC_VCF=$(echo $RM_DUP_VCF | sed 's/.vcf/_fixgt.vcf.gz/') # append _fixtgt
    local BEST_IMPUTE_VCF=$(echo $VCF | sed 's/.vcf.gz/_impute/') # append _impute and remove .vcf.gz 
    # PreProcessing
    #
    singularity exec ${SIF_DIR}/htslib%3A1.18--h81da01d_0 tabix -f -p vcf $VCF # (re)write .tbi
    singularity exec ${SIF_DIR}/bcftools%3A1.17--haef29d1_0 bcftools norm -f $GENOME -Ov -o ${RM_DUP_VCF} -d none $VCF # remove duplicates
    singularity exec ${SIF_DIR}/bcftools%3A1.17--haef29d1_0 sh -c "bcftools +fixploidy $RM_DUP_VCF | bcftools view -Oz -o $BIALLELIC_VCF" # fix ploidy
    singularity exec ${SIF_DIR}/htslib%3A1.18--h81da01d_0 tabix -f -p vcf $BIALLELIC_VCF # write .tbi
    singularity exec ${SIF_DIR}/beagle%3A4.1_21Jan17.6cc.jar--0 beagle -Xmx10g gt=$BIALLELIC_VCF out=$BEST_IMPUTE_VCF impute=false # phase genotypes
    singularity exec ${SIF_DIR}/htslib%3A1.18--h81da01d_0 tabix -f -p vcf "${BEST_IMPUTE_VCF}.vcf.gz" # write .tbi
    echo "${BEST_IMPUTE_VCF}.vcf.gz"
}
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

    # Run sequenceDiversity for each population individually
    for i in $(seq 1 $NUM_POPULATIONS); do
        if [ ! -f "$OUTDIR/${prefix}_pop${i}_pi.txt" ]; then
            singularity exec ${SIF_DIR}/vcflib%3A1.0.3--hecb563c_1 sequenceDiversity --target "${pop_indices_array[$i]}" -y GT --file "${VCF}" > "$OUTDIR/gwas/${prefix}_pop${i}_pi.txt"
        fi
    done
}

wait_for_job_to_finish() {
    JOBID=$1
    while squeue -j $JOBID | grep -q $JOBID; do
        sleep 10  # wait for 10 seconds before checking again
    done
    echo "Job $JOBID has finished."
}

full_simulations_workflows() {
    local OUTDIR="$1"
    local popsize="${2:-20000}" # second argument defaults to 20000
    local n_pop="${3:-4}" # third argument defaults to 3
    local n_indv_per_pop="${4:-30}" # fourth argument defaults to 4
    local CHROMOSOME_LIST="$5"
    local DENOVO="${6:-'false'}"
    local DIRNAME=$(basename $OUTDIR) # PREFIX
    local TRUTH_VCF="${OUTDIR}/ref_loci_vars/ri_master.vcf.gz" # FILE
    local TRUTH_FASTA="${OUTDIR}/ref_loci_vars/reference_rad_loci.fa.gz" # FILE
    local TRUTH_BGZIP_VCF=$(echo $TRUTH_VCF | sed 's/ri_master/ri-master_bgzip/') # FILE
    local TRUTH_SUB_VCF=$(echo $TRUTH_VCF | sed 's/ri_master/ri-master_subset/') # FILE
    local TRUTH_SUB_NORM_VCF=$(echo $TRUTH_SUB_VCF | sed 's/subset/subset_norm/') # FILE
    local TRUTH_REGIONS=$(echo $TRUTH_VCF | sed 's/.vcf.gz/.bed/') # FILE
    local RAD_READS_DIR=$(echo $OUTDIR | sed 's,$,/rad_reads,') # DIRECTORY
    # REFERENCE LINEAR
    local NF_RADSEQ_REFERENCE_BED="${OUTDIR}/nf-radseq/reference/bwa-mem2/intervals/bedops_merge/*.bed" # FILE
    local RADSEQ_REF_VCF_DIR="${OUTDIR}/nf-radseq/reference/variant_calling" # DIRECTORY
    local RAD_REF_VCF_DIR="${RADSEQ_REF_VCF_DIR}/filter" # DIRECTORY
    local RAD_REF_VCF_NORM_FILE="${RADSEQ_REF_VCF_DIR}/msp_norm.vcf.gz"
    local RAD_REF_VCF_NORM_FILTERED_FILE="${RADSEQ_REF_VCF_DIR}/msp_norm_qual20.vcf.gz"
    local RAD_REF_VCF_NORM_FILTERED_MAF_FILE="$(echo $RAD_REF_VCF_NORM_FILTERED_FILE | sed 's/.vcf.gz/_maf01.vcf.gz/g')"
    # UNIVERSAL VARIANT GRAPH VARIABLES
    local FQ="${RAD_READS_DIR}/*.{1,2}.fq.gz"
    # REFERENCE VARIANT GRAPH
    local REF_FAI="${RADSEQ_REF_VCF_DIR}/../samtools/index/*.fasta.fai"
    local RAD_VG_REF_VCF_DIR="${OUTDIR}/vg/reference/variant_calling/filter"
    local RAD_VG_REF_VCF_NORM_FILE="${OUTDIR}/vg/reference/BCFTOOLS/NORM/msp_norm.vcf.gz"
    local RAD_VG_REF_VCF_NORM_FILTERED_FILE="${OUTDIR}/vg/reference/BCFTOOLS/NORM/msp_norm_qual20.vcf.gz"
    local RAD_VG_REF_VCF_NORM_FILTERED_FILE_AF=$(echo "$RAD_VG_REF_VCF_NORM_FILTERED_FILE" | sed 's/.vcf.gz/_af.vcf.gz/')
    # DE NOVO LINEAR
    local RADSEQ_DENOVO_VCF_DIR="${OUTDIR}/nf-radseq/denovo/variant_calling" # DIRECTORY
    local NF_RADSEQ_DENOVO_BED=$(echo $NF_RADSEQ_REFERENCE_BED | sed 's/reference/denovo/') # FILE
    local RAD_DENOVO_VCF_DIR=$(echo $RAD_REF_VCF_DIR | sed 's/reference/denovo/')
    local RAD_VG_DENOVO_VCF_DIR=$(echo $RAD_VG_REF_VCF_DIR | sed 's/reference/denovo/')

    #
    # SIMULATIONS
    #
    echo " Parameters:"
    echo "  Reference sequences: ${GENOME}"
    echo "  Chromosome list: ${CHROMOSOME_LIST}"
    echo "  truth vcf: $TRUTH_VCF"

    if [ ! -d "$OUTDIR" ] || [ -z "$(ls -A "$OUTDIR")" ]; then
        mkdir -p "$OUTDIR/test_accuracy"
        echo "Writing radinitio output to: $OUTDIR"
        # Simulating a ddRAD library:
        module -q load python/3.9.1
        radinitio --simulate-all \
            --genome $GENOME \
            --chromosomes $CHROMOSOME_LIST \
            --out-dir $OUTDIR \
            --n-pops $n_pop --pop-eff-size $popsize --n-seq-indv $n_indv_per_pop \
            --library-type ddRAD --enz PstI --enz2 EcoRI \
            --insert-mean 400 \
            --pcr-cycles 5 --coverage 20 --read-length 150
    else
        echo "Directory already exists: $OUTDIR"
    fi

    singularity exec ${SIF_DIR}/htslib%3A1.18--h81da01d_0 sh -c "bgzip -d -c $TRUTH_FASTA | grep '>' | cut -f 2 -d'=' | sed 's/[:-]/\t/g' | mawk '{x=$2+5;y=$3-5; print $1"\t"x"\t"y}' > ${TRUTH_REGIONS}" # create bed file
    singularity exec ${SIF_DIR}/bcftools%3A1.17--haef29d1_0 bcftools view -Oz -o $TRUTH_BGZIP_VCF $TRUTH_VCF
    singularity exec ${SIF_DIR}/htslib%3A1.18--h81da01d_0 tabix -p vcf -f $TRUTH_BGZIP_VCF # write .tbi
    singularity exec ${SIF_DIR}/bcftools%3A1.17--haef29d1_0 bcftools view -R ${TRUTH_REGIONS} -Ov  $TRUTH_BGZIP_VCF | vcftools --vcf - --mac 3 --recode --recode-INFO-all --stdout| bcftools view -Ov -o $TRUTH_SUB_VCF #You need to match the default filterinf in RADSeq 
    singularity exec ${SIF_DIR}/htslib%3A1.18--h81da01d_0 tabix -p vcf -f $TRUTH_SUB_VCF # write .tbi
    singularity exec ${SIF_DIR}/bcftools%3A1.17--haef29d1_0 bcftools norm -a -m -any -f $GENOME $TRUTH_SUB_VCF > $TRUTH_SUB_NORM_VCF
    singularity exec ${SIF_DIR}/bcftools%3A1.17--haef29d1_0 bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $TRUTH_SUB_NORM_VCF > "${OUTDIR}/test_accuracy/ri-master_norm_dedup_rename_query.txt"

    extract_population_indices_and_calculate "$TRUTH_SUB_VCF" "$OUTDIR" 'ri-master' "$n_pop" # actual pop. stats
    # Switch python versions
    module -q load python/3.7.13 # required in order for hap.py to work
    #
    # Write dummy .fq.gz files from .fa.gz radinition output
    #
    echo "writing fastq files with dummy quality scores (@)"
    for fa in "${RAD_READS_DIR}"/*.fa.gz; do
        fq=$(echo $fa | sed 's/.fa.gz/.fq.gz/')
        if [ -f "$fq" ] && [ -s "$fq" ]; then
            echo "$fq exists and is not empty."
        else
            #echo "$fq does not exist or is empty."
            singularity exec ${SIF_DIR}/seqtk%3A1.4--he4a0461_1 seqtk seq -F '@' $fa | gzip -c > $fq &
            # Optional: limit the number of background jobs
            if (( i % 15 == 0 )); then
                wait
            fi
        fi
    done
    wait
    echo "finished writing fastq files with dummy quality scores (@)"
    #
    # Write Samplesheet for nf-core/radseq
    #
    if [ ! -f ${RAD_READS_DIR}/$DIRNAME.csv ]; then 
        echo "Writing samplesheet: ${RAD_READS_DIR}/$DIRNAME.csv"
        # Write SampleSheet
        echo "sample,fastq_1,fastq_2,umi_barcodes" > "${RAD_READS_DIR}/$DIRNAME.csv"
        paste -d',' <(for i in "${RAD_READS_DIR}"/*.1.fq.gz; do basename $i | cut -f1 -d'.' -; done)\
        <(ls "${RAD_READS_DIR}"/*.1.fq.gz) <(ls "${RAD_READS_DIR}"/*.2.fq.gz)\
        <(for i in "${RAD_READS_DIR}"/*.1.fq.gz; do echo "false"; done)\
        >> "${RAD_READS_DIR}/$DIRNAME.csv"
    else
        echo "samplesheet: ${RAD_READS_DIR}/$DIRNAME.csv already exists. Skipping."
    fi
    SAMPLESHEET="${RAD_READS_DIR}/$DIRNAME.csv"
    
    #
    # Run NextFlow radseq reference
    #
    JOBID=$(sbatch --job-name 'reference_radseq' --time=72:00:00 -p uri-cpu --mem=2G -c 1 -e %x_%j.err -o %x_%j.out ../../../../nextflow_workflows_stache.sh "$OUTDIR" 'radseq' 'reference' "$SAMPLESHEET" "$GENOME" | awk '{print $4}')
    #bash ../../../nextflow_workflows_stache.sh "$OUTDIR" 'radseq' 'reference' "$SAMPLESHEET" "$GENOME"
    wait_for_job_to_finish $JOBID
    # Write chromosome and position file across reference vcfs
    echo "creating _chrompos.txt inside ${RAD_REF_VCF_DIR}"
    for vcf in ${RAD_REF_VCF_DIR}/*.vcf.gz; do
        name=$(basename $vcf | cut -f 1,2 -d '.')
        singularity exec ${SIF_DIR}/bcftools%3A1.17--haef29d1_0 bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $vcf > "${RAD_REF_VCF_DIR}/${name}_chrompos.txt" &
    done
    wait
    #BEST_REF_VCF=$(paste -d '' <(echo "${RAD_REF_VCF_DIR}/") <(python3 ../../../../findMostAccurateVCF.py "${RAD_REF_VCF_DIR}" True)) # grab linear-based variant file w/ highest f1 score
    singularity exec ${SIF_DIR}/bcftools%3A1.17--haef29d1_0 bcftools view -i 'QUAL>20' -Oz -o $RAD_REF_VCF_NORM_FILTERED_FILE $RAD_REF_VCF_NORM_FILE
    BEST_REF_PHASED_VCF=$(vcf_rmdup_fixgt_phase $RAD_REF_VCF_NORM_FILTERED_FILE | tail -n 1)
    singularity exec ${SIF_DIR}/hap.py%3A0.3.15--py27hcb73b3d_0 hap.py "$TRUTH_SUB_NORM_VCF" "$BEST_REF_PHASED_VCF" -r $GENOME -o "${OUTDIR}/test_accuracy/nf-radseq-reference"
    
    extract_population_indices_and_calculate $BEST_REF_PHASED_VCF $OUTDIR "nf-radseq-reference" 
    singularity exec ${SIF_DIR}/bcftools%3A1.17--haef29d1_0 bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' $RAD_REF_VCF_NORM_FILTERED_FILE > "${OUTDIR}/test_accuracy/nf-radseq-reference_chrompos.txt"
    singularity exec ${SIF_DIR}/bcftools%3A1.17--haef29d1_0 bcftools view -i 'MAF>0.01' -Oz -o "$RAD_REF_VCF_NORM_FILTERED_MAF_FILE" "$RAD_REF_VCF_NORM_FILTERED_FILE" 

    singularity exec ${SIF_DIR}/htslib%3A1.18--h81da01d_0 tabix -p vcf -f $RAD_REF_VCF_NORM_FILTERED_MAF_FILE # radseq doesn't index reference
    REF_TBI="${RAD_REF_VCF_NORM_FILTERED_MAF_FILE}.tbi"
    echo $RAD_REF_VCF_NORM_FILTERED_MAF_FILE
    #
    # Variant Graph
    #
    VG_JOBID=$(sbatch --dependency=afterok:$JOBID --job-name 'reference_vg' --time=72:00:00 -p uri-cpu --mem=5G -c 1 -e %x_%j.err -o %x_%j.out ../../../../nextflow_workflows_stache.sh "$OUTDIR/vg/reference" 'vg' "$GENOME" "$REF_FAI" "$RAD_REF_VCF_NORM_FILTERED_MAF_FILE" "$REF_TBI" "$FQ" | awk '{print $4}')
    wait_for_job_to_finish $VG_JOBID
    for vcf in ${RAD_VG_REF_VCF_DIR}/*.vcf.gz; do
        name=$(basename $vcf | cut -f 1,2 -d '.')
        singularity exec ${SIF_DIR}/bcftools%3A1.17--haef29d1_0 bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $vcf > "${RAD_VG_REF_VCF_DIR}/${name}_chrompos.txt" &
    done
    wait
    #BEST_VG_REF_VCF=$(paste -d '' <(echo "${RAD_VG_REF_VCF_DIR}/") <(python3 ../../../../findMostAccurateVCF.py "${RAD_VG_REF_VCF_DIR}" True))
    singularity exec ${SIF_DIR}/bcftools%3A1.17--haef29d1_0 bcftools view -i 'QUAL>20' -Oz -o $RAD_VG_REF_VCF_NORM_FILTERED_FILE $RAD_VG_REF_VCF_NORM_FILE
    singularity exec ${SIF_DIR}/bcftools%3A1.17--haef29d1_0 bcftools +fill-tags "$RAD_VG_REF_VCF_NORM_FILTERED_FILE" -Oz -o "$RAD_VG_REF_VCF_NORM_FILTERED_FILE_AF" -- -t 'AF'
    singularity exec ${SIF_DIR}/bcftools%3A1.17--haef29d1_0 bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$RAD_VG_REF_VCF_NORM_FILTERED_FILE_AF" > "${OUTDIR}/test_accuracy/nf-vg-reference_chrompos.txt"
    singularity exec "${SIF_DIR}/htslib%3A1.18--h81da01d_0" tabix -f -p vcf "${RAD_VG_REF_VCF_NORM_FILTERED_FILE}"
    BEST_VG_REF_PHASED_VCF=$(vcf_rmdup_fixgt_phase $RAD_VG_REF_VCF_NORM_FILTERED_FILE | tail -n 1)
    singularity exec ${SIF_DIR}/hap.py%3A0.3.15--py27hcb73b3d_0 hap.py "$TRUTH_SUB_NORM_VCF" "$BEST_VG_REF_PHASED_VCF" -r $GENOME -o $OUTDIR/test_accuracy/nf-vg-reference
    extract_population_indices_and_calculate $BEST_VG_REF_PHASED_VCF $OUTDIR "nf-vg-reference" 

    # move files into test_accuracy
    #BEST_REF_VCF_CHROMPOS=$(echo "${RAD_REF_VCF_NORM_FILTERED_FILE_CHROMPOS}" | sed 's/.vcf.gz/_chrompos.txt/')
    cp "${RAD_REF_VCF_DIR}/msp_norm_qual20_chrompos.txt" "${OUTDIR}/test_accuracy/nf-radseq-reference_chrompos.txt"
    BEST_VG_REF_VCF_CHROMPOS=$(echo "${BEST_VG_REF_VCF}" | sed 's/.vcf.gz/_chrompos.txt/')
    grep -v '*' "${RAD_VG_REF_VCF_DIR}/msp_norm_qual20_chrompos.txt" > "${OUTDIR}/test_accuracy/nf-vg-reference_chrompos.txt"
    #cp ${BEST_VG_REF_VCF_CHROMPOS} "${OUTDIR}/test_accuracy/nf-vg-reference_chrompos.txt"

    if [[ $DENOVO == true ]]; then
        #
        # Run radseq
        #
        DENOVO_RADSEQ_JOBID=$(sbatch --job-name 'denovo_radseq' --time=72:00:00 -p uri-cpu --mem=2G -c 1 -e %x_%j.err -o %x_%j.out ../../../../nextflow_workflows_stache.sh "$OUTDIR" 'radseq' 'denovo' "$SAMPLESHEET" '2,3' '2,3' | awk '{print $4}')
        wait_for_job_to_finish $DENOVO_RADSEQ_JOBID
        # Write chromosome and position file across reference vcfs
        echo "creating _chrompos.txt inside $RAD_DENOVO_VCF_DIR"
        for vcf in ${RAD_DENOVO_VCF_DIR}/*.vcf.gz; do
            name=$(basename $vcf | cut -f 1,2 -d '.')
            singularity exec "${SIF_DIR}/bcftools%3A1.17--haef29d1_0" bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $vcf > "${RAD_DENOVO_VCF_DIR}/${name}_chrompos.txt" &
        done
        wait
        BEST_DENOVO_VCF=$(paste -d '' <(echo "${RAD_DENOVO_VCF_DIR}/") <(python3 ../../../../findMostAccurateVCF.py "${RAD_DENOVO_VCF_DIR}/" True))
        BEST_DENOVO_PHASED_VCF=$(vcf_rmdup_fixgt_phase $BEST_DENOVO_VCF | tail -n 1)
        extract_population_indices_and_calculate $BEST_DENOVO_PHASED_VCF $OUTDIR "nf-radseq-denovo"
        #
        # Define variables for vg in denovo mode
        #
        DENOVO_FASTA=$(echo $BEST_DENOVO_VCF | sed -E 's/msp_([0-9]+)_([0-9]+)_.+\.vcf\.gz/msp_\1_\2_rainbow.fasta/' | sed 's,variant_calling/filter/,write_fasta/,') # take the denovo vcf outputted by BEST_DENOVO_VCF and use denovo combination to select reference
        DENOVO_FAI=$(paste -d '' <(echo "$OUTDIR/nf-radseq/denovo/samtools/index/") <(basename $BEST_DENOVO_VCF | sed -E 's/msp_([0-9]+)_([0-9]+)_.+\.vcf\.gz/msp_\1_\2_rainbow.fasta/') <(echo '.fai'))
        singularity exec "${SIF_DIR}/htslib%3A1.18--h81da01d_0" tabix -f -p vcf ${BEST_DENOVO_VCF} # radseq doesn't index reference
        DENOVO_TBI="${BEST_DENOVO_VCF}.tbi"
        #
        # Run vg
        #
        DENOVO_VG_JOBID=$(sbatch --dependency=afterok:$DENOVO_RADSEQ_JOBID --job-name 'denovo_vg' --time=72:00:00 -p uri-cpu --mem=2G -c 1 -e %x_%j.err -o %x_%j.out ../../../../nextflow_workflows_stache.sh "$OUTDIR/vg/denovo" 'vg' "$DENOVO_FASTA" "$DENOVO_FAI" "$BEST_DENOVO_VCF" "$DENOVO_TBI" "$FQ" | awk '{print $4}')
        wait_for_job_to_finish $DENOVO_VG_JOBID
        RAD_VG_DENOVO_VCF_DIR=$(echo $RAD_VG_REF_VCF_DIR | sed 's/reference/denovo/')
        for vcf in ${RAD_VG_DENOVO_VCF_DIR}/*.vcf.gz; do
            name=$(basename $vcf | cut -f 1,2 -d '.')
            singularity exec "${SIF_DIR}/bcftools%3A1.17--haef29d1_0" bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $vcf > "${RAD_VG_DENOVO_VCF_DIR}/${name}_chrompos.txt"
        done
        BEST_VG_DENOVO_VCF=$(paste -d '' <(echo "${RAD_VG_DENOVO_VCF_DIR}/") <(python3 ../../../../findMostAccurateVCF.py "${RAD_VG_DENOVO_VCF_DIR}" True))
        echo "${BEST_VG_DENOVO_VCF}"
        singularity exec "${SIF_DIR}/htslib%3A1.18--h81da01d_0" -f -p vcf "${BEST_VG_DENOVO_VCF}"

        BEST_VG_DENOVO_PHASED_VCF=$(vcf_rmdup_fixgt_phase $BEST_VG_DENOVO_VCF | tail -n 1)
        extract_population_indices_and_calculate $BEST_VG_DENOVO_PHASED_VCF $OUTDIR "nf-vg-denovo"
    fi
}
#######################################################################################################################################################
#                                                                MAIN WORKFLOW
# ----------------------
# INPUT ARGUMENTS
# ----------------------
# Set Defaults
GENOME="$(pwd)/atlantic_hsc/Lp_genome_Chr26.fasta"
n_pop=4
n_indv_per_pop=30
type="low_gene_flow"
outdir=$(pwd) # Default output directory is the current working directory

# Parse User Parameters
while [ "$1" != "" ]; do
    case $1 in
        -g | --genome )                 shift
                                        GENOME=$1
                                        ;;
        -n | --number_of_populations )  shift
                                        n_pop=$1
                                        ;;
        -i | --number_of_individuals )  shift
                                        n_indv_per_pop=$1
                                        ;;
        -o | --outdir )                 shift
                                        outdir=$1
                                        ;;
        -t | --type )                   shift
                                        type=$1
                                        ;;
        * )                     echo "Invalid parameter detected: $1"
                                exit 1
                                ;;
    esac
    shift
done

# After parsing parameters, you can set CHROMOSOME_LIST
CHROMOSOME_LIST=$(echo $GENOME | sed 's/\(\.fa\|\.fasta\)$/_chromlist.txt/') # file matching sequence headers in fasta to do simulations

for effpopsize in 5000 15000; do
    for replicate in $(seq 1 30); do
        
        # Simulation Variables
        OUTDIR=$(echo $GENOME | sed "s,\(\.fa\|\.fasta\)$,_${type}/${effpopsize}/${replicate}/,") # DIRECTORY
        mkdir -p $OUTDIR
        cd $OUTDIR
        
        full_simulations_workflows "$OUTDIR" "$effpopsize" "$n_pop" "$n_indv_per_pop" "$CHROMOSOME_LIST" false &
        cd ../../../ # back out of ./atlantic_hsc/Lp_genome_Chr26_10000/1
        if (( ++count % 6 == 0 )); then # limit the number of simulation workflows running at once
            wait
        fi
    done
    wait
done


