#!/bin/bash
export NXF_WORK="/scratch/workspace/gabriel_barrett_uri_edu-work/sims"
 OUTDIR="$1"
# Check the first input argument
if [[ "$2" == "vg" ]]; then
    echo "You chose vg!"
    if [[ $# -ne 7 ]]; then
        echo "Error: When choosing 'vg', you must provide 5 arguments in total.
                radseq_workflows.sh vg <.fasta> <.fai> <.vcf> <.tbi> <.fq.gz>"
        exit 1
    else
        echo "You chose vg with correct number of arguments!"
         REFERENCE_FASTA="$3"
         REFERENCE_FAI="$4"
         REFERENCE_VCF="$5"
         REFERNCE_TBI="$6"
         FQ="$7"
        
        nextflow run Gabriel-A-Barrett/nf-vg-pipeline --fasta ${REFERENCE_FASTA} --fai "${REFERENCE_FAI}" \
        --vcf "${REFERENCE_VCF}" --tbi "${REFERNCE_TBI}" -resume \
        --fq "${FQ}" --output_mode 'symlink' -c ~/config/unity.config --outdir "${OUTDIR}" --remove_duplicates 'true' --atomize_variants 'true'

    fi

elif [[ "$2" == "radseq" ]]; then
    MODE="$3"
    SAMPLESHEET="$4"
    if [[ "$MODE" == "reference" ]]; then
        echo "You chose radseq reference"
        if [[ $# -ne 5 ]]; then
        echo "Error: When choosing 'radseq' under reference mode, you must provide 4 arguments in total.
                radseq_workflows.sh 'radseq' <MODE> <SAMPLESHEET> <REFERENCE_FASTA>"
        exit 1
        fi
        REFERENCE_FASTA="$5"
        nextflow run Gabriel-A-Barrett/radseq -r dev --input "$SAMPLESHEET" -c ~/config/unity.config --sequence_type 'PE'\
         --method 'reference' --genome "$REFERENCE_FASTA" --outdir "${OUTDIR}/nf-radseq" -resume --save_reference_indices 'true' --atomize_variants 'true' --save_trimmed='true' --use_best_n_alleles='false'
            
    elif [[ "$MODE" == "denovo" ]]; then
        echo "Using denovo mode for radseq."
        MIN_READS_WITHIN_INDV="${5:-2,3}"
        MIN_READS_BETWEEN_INDV="${6:-2,3}"
        nextflow run Gabriel-A-Barrett/radseq -r dev --input "$SAMPLESHEET" -c ~/config/unity.config --sequence_type 'PE'\
        --method 'denovo'  --outdir "${OUTDIR}/nf-radseq" -resume --save_reference_indices 'true' --atomize_variants 'true'\
        --minreaddepth_withinindividual="$MIN_READS_WITHIN_INDV" --minreaddepth_betweenindividual="$MIN_READS_BETWEEN_INDV" --only_use_best_reference 'true'
            
    else
        echo "Error: Mode should be either 'reference' or 'denovo' for radseq."
        exit 1
    fi
else
    echo "Invalid input. Please choose either 'vg' or 'radseq'."
fi
