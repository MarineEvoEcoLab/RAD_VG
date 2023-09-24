#!/usr/bin/env bash

set -e # Exit on error

module load apptainer/1.1.5+py3.8.12
#module load python/3.9.1
#module load bcftools/1.15
#module load htslib/1.14
#module load parallel/20210922

export NXF_SINGULARITY_CACHEDIR="./nxf-singularity-cache-dir"

calcTruthPopStats() {
    local TRUTH_NORM_VCF=$1
    prefix=$(echo $TRUTH_NORM_VCF | sed 's/.vcf.gz//')
    local OUTDIR=$2

    ref_pop0=$(seq -s, 1 10)
    ref_pop1=$(seq -s, 11 20)
    sequence_diversity_pop=$(seq -s, 1 20)
    singularity exec vcflib%3A1.0.3--hecb563c_1 pFst --target $ref_pop0 --background $ref_pop1 -y GT --file $TRUTH_NORM_VCF > ${OUTDIR}/ri-master_pFst.txt
    singularity exec vcflib%3A1.0.3--hecb563c_1 sequenceDiversity --target $sequence_diversity_pop -y GT --file $TRUTH_NORM_VCF > ${OUTDIR}/ri-master_pi.txt

}
writeTruthDataFiles() {

    local OUTDIR="${1}"
    local TRUTH_VCF="${OUTDIR}/ref_loci_vars/ri-master.vcf.gz" # FILE
    local TRUTH_NORM_VCF=$(echo $TRUTH_VCF | sed 's/ri-master/ri-master_norm_rename_subset/') # FILE

    gunzip -c ${OUTDIR}/ref_loci_vars/reference_rad_loci.fa.gz | grep '>' | cut -f 2 -d'=' | sed 's/[:-]/\t/g' > ${OUTDIR}/ref_loci_vars/truth.bed

    singularity exec bcftools%3A1.17--haef29d1_0 bcftools query -l $TRUTH_VCF | awk -v OFS='\t' '{print $1,$2=$1"_ref"}' > ${OUTDIR}/ref_loci_vars/rename_reference.txt
    
    gunzip $TRUTH_VCF 
    bgzip $(echo $TRUTH_VCF | sed 's/.gz//')  

    singularity exec bcftools%3A1.17--haef29d1_0 bcftools norm -d all -f $GENOME $TRUTH_VCF | bcftools reheader -s ${OUTDIR}/ref_loci_vars/rename_reference.txt > ${OUTDIR}/ref_loci_vars/ri-master_norm_rename.vcf
    singularity exec htslib%3A1.18--h81da01d_0 bgzip ${OUTDIR}/ref_loci_vars/ri-master_norm_rename.vcf
    singularity exec htslib%3A1.18--h81da01d_0 tabix -f -p vcf ${OUTDIR}/ref_loci_vars/ri-master_norm_rename.vcf.gz
    singularity exec bcftools%3A1.17--haef29d1_0 bcftools view -R ${OUTDIR}/ref_loci_vars/truth.bed -Oz -o $TRUTH_NORM_VCF ${OUTDIR}/ref_loci_vars/ri-master_norm_rename.vcf.gz
    singularity exec bcftools%3A1.17--haef29d1_0 bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $TRUTH_NORM_VCF > ${OUTDIR}/ref_loci_vars/../test_accuracy/ri-master_norm_dedup_rename_query.txt 
    singularity exec htslib%3A1.18--h81da01d_0 tabix -f -p vcf $TRUTH_NORM_VCF

    calcTruthPopStats $TRUTH_NORM_VCF $OUTDIR
}
calculate_population_statistics() {
    # Calculate Population Statistics For reference radseq
    if [ -z "$1" ]; then
        echo "Usage: calculate_population_statistics <VCF>"
        return 1
    fi

    VCF="$1"
    TRUTH_NORM_VCF="$2"
    prefix="$3"
    OUTDIR="$4"
    RM_DUP_VCF=$(echo $VCF | sed 's/.vcf.gz/_dedup.vcf/')
    BIALLELIC_VCF=$(echo $RM_DUP_VCF | sed 's/.vcf/_fixgt.vcf.gz/')
    BED=$(echo $VCF | sed 's/.vcf.gz/.bed/')
    TRUTH_SUBSET_VCF=$(echo $TRUTH_NORM_VCF | sed "s/.vcf.gz/_${prefix}.vcf.gz/")
    TEST_SUBSET_VCF=$(echo $BIALLELIC_VCF | sed "s/.vcf/_subset.vcf/")

    if [[ $prefix == *denovo* ]]; then
        echo "calculating population statistics based on model vcf"
        pop0_indvs=$(awk 'NR <= 10 { next } { printf("%s%s", sep, $1); sep="," } END { printf("\n") }' $BEST_IMPUTE_VCF_INDVS)
        pop1_indvs=$(awk 'NR > 10 { next } { printf("%s%s", sep, $1); sep="," } END { printf("\n") }' $BEST_IMPUTE_VCF_INDVS)
        sequence_diversity_pop=$(seq -s, 1 20)
        singularity exec vcflib%3A1.0.3--hecb563c_1 pFst --target $pop0_indvs --background $pop1_indvs -y GT --file "${BEST_IMPUTE_VCF}.vcf.gz" > $OUTDIR/${prefix}_pFst.txt
        singularity exec vcflib%3A1.0.3--hecb563c_1 sequenceDiversity --target $sequence_diversity_pop -y GT --file "${BEST_IMPUTE_VCF}.vcf.gz" > $OUTDIR/${prefix}_pi.txt

    else
        singularity exec htslib%3A1.18--h81da01d_0 tabix -f -p vcf $VCF
        singularity exec bcftools%3A1.17--haef29d1_0 sh -c "bcftools norm -Ov -o ${RM_DUP_VCF} -d none $VCF && tabix -p vcf ${RM_DUP_VCF}.gz"
        echo "fixing ploidy"
        singularity exec bcftools%3A1.17--haef29d1_0 sh -c "bcftools +fixploidy $RM_DUP_VCF | bcftools view -Oz -o $BIALLELIC_VCF && tabix -f -p vcf $BIALLELIC_VCF"
        echo "converting vcf to bed file"
        # Converting $UNZIP_REF_VCF to $BED
        singularity exec bedops%3A2.4.41--h9f5acd7_0 vcf2bed --max-mem=2G < $RM_DUP_VCF > $BED
        echo "created regions bed file: $BED"
        # Subset reference to regions
        singularity exec htslib%3A1.18--h81da01d_0 tabix -f -p vcf $TRUTH_NORM_VCF
        # Get reference snps for sites in model
        singularity exec bcftools%3A1.17--haef29d1_0 bcftools view -R $BED -Oz -o $TRUTH_SUBSET_VCF $TRUTH_NORM_VCF
        #singularity exec htslib%3A1.18--h81da01d_0 tabix -f -p vcf $TRUTH_SUBSET_VCF        
        # Now reduce vg vcf down with
        singularity exec bcftools%3A1.17--haef29d1_0 bcftools isec -n +2 -w 1 -o $TEST_SUBSET_VCF $TRUTH_SUBSET_VCF $BIALLELIC_VCF
        # Impute vcf for genotype phasing for vcflib and write genotype probabilities
        BEST_IMPUTE_VCF=$(echo $VCF | sed 's/.vcf.gz/_impute/')
        echo "writing imputed genotypes to ${BEST_IMPUTE_VCF}.vcf.gz"
        
        #
        # Imputation w/ Beagle
        #
        singularity exec beagle%3A5.2_21Apr21.304--hdfd78af_0 beagle -Xmx10g gt=$BIALLELIC_VCF out=$BEST_IMPUTE_VCF ref=$TRUTH_SUBSET_VCF gp=true
        singularity exec htslib%3A1.18--h81da01d_0 tabix -f -p vcf "${BEST_IMPUTE_VCF}.vcf.gz"

        # Calculate window-based pi and fst
        # How to get the sample order
        BEST_IMPUTE_VCF_INDVS=$(echo "${BEST_IMPUTE_VCF}.vcf.gz" | sed 's/.vcf.gz/_indv_order.tsv/')
        singularity exec bcftools%3A1.17--haef29d1_0 bcftools query -l "${BEST_IMPUTE_VCF}.vcf.gz" | awk -v OFS='\t' '{print NR, $0}' | sort -k2 > $BEST_IMPUTE_VCF_INDVS
        echo $BEST_IMPUTE_VCF_INDVS
        pop0_indvs=$(awk 'NR <= 10 { next } { printf("%s%s", sep, $1); sep="," } END { printf("\n") }' $BEST_IMPUTE_VCF_INDVS)
        pop1_indvs=$(awk 'NR > 10 { next } { printf("%s%s", sep, $1); sep="," } END { printf("\n") }' $BEST_IMPUTE_VCF_INDVS)
        sequence_diversity_pop=$(seq -s, 1 20)
        
        echo $pop0_indvs
        echo $pop_indvs

        singularity exec vcflib%3A1.0.3--hecb563c_1 pFst --target $pop0_indvs --background $pop1_indvs -y GT --file "${BEST_IMPUTE_VCF}.vcf.gz" > $OUTDIR/${prefix}_pFst.txt
        singularity exec vcflib%3A1.0.3--hecb563c_1 sequenceDiversity --target $sequence_diversity_pop -y GT --file "${BEST_IMPUTE_VCF}.vcf.gz" > $OUTDIR/${prefix}_pi.txt
        
        ref_pop0=$(seq -s, 1 10)
        ref_pop1=$(seq -s, 11 20)
        singularity exec htslib%3A1.18--h81da01d_0 tabix -f -p vcf ${TRUTH_SUBSET_VCF}
        singularity exec vcflib%3A1.0.3--hecb563c_1 pFst --target $ref_pop0 --background $ref_pop1 -y GT --file "${TRUTH_SUBSET_VCF}" > $OUTDIR/${prefix}_truth_pFst_tmp.txt
        singularity exec vcflib%3A1.0.3--hecb563c_1 sequenceDiversity --target $sequence_diversity_pop -y GT --file "${TRUTH_SUBSET_VCF}" > $OUTDIR/${prefix}_truth_pi.txt

        # adding removal of extra variants from truth dataset
        awk '{print $1,$2}' $OUTDIR/${prefix}_pFst.txt > $OUTDIR/${prefix}_target_pFst.txt
        awk 'NR==FNR {a[$1,$2]++; next} (a[$1,$2])' $OUTDIR/${prefix}_target_pFst.txt $OUTDIR/${prefix}_truth_pFst_tmp.txt > $OUTDIR/${prefix}_truth_pFst.txt
        rm $OUTDIR/${prefix}_target_pFst.txt $OUTDIR/${prefix}_truth_pFst_tmp.txt
        echo "Created genome-wide statistics for: $OUTDIR/${prefix}"
    fi
}
split_variant_classes() {
    declare VCF="$1" OUTDIR="$2"
    prefix=$(basename $VCF | sed 's/.vcf.gz//')    
    SNP_OUTPUT="${OUTDIR}/test_accuracy/${prefix}_snps.vcf.gz"
    SV_OUTPUT="${OUTDIR}/test_accuracy/${prefix}_sv.vcf.gz"
    echo "splitting VCF based on variant class"
    singularity exec bcftools%3A1.17--haef29d1_0 bcftools view -Oz -o "$SNP_OUTPUT" -v snps,mnps $VCF &
    singularity exec bcftools%3A1.17--haef29d1_0 bcftools view -Oz -o "$SV_OUTPUT" -v indels $VCF &  
    wait
    echo "split_variant_classes: ${VCF} processing completed"
}

GENOME=$1
echo "Running Simulations based on sequence(s): $1"
popsize="${2:-20000}" # second argument defaults to 20000
CHROMOSOME_LIST=$(echo $GENOME | sed 's/\(\.fa\|\.fasta\)$/_chromlist.txt/')
OUTDIR=$(echo $GENOME | sed "s/\(\.fa\|\.fasta\)$/_$popsize/") # DIRECTORY
DIRNAME=$(basename $OUTDIR) # PREFIX
TRUTH_VCF="${OUTDIR}/ref_loci_vars/ri-master.vcf.gz" # FILE
TRUTH_NORM_VCF=$(echo $TRUTH_VCF | sed 's/ri-master/ri-master_norm_rename_subset/') # FILE

export NXF_TEMP="${OUTDIR}/work"
export NXF_WORK="${OUTDIR}/work"

if [ ! -d "$OUTDIR" ] || [ -z "$(ls -A "$OUTDIR")" ]; then
    mkdir -p "$OUTDIR/test_accuracy"
    echo "Writing files to: $OUTDIR"
    # Simulating a ddRAD library:
    radinitio --simulate-all \
        --genome $GENOME \
        --chromosomes $CHROMOSOME_LIST \
        --out-dir $OUTDIR \
        --n-pops 2 --pop-eff-size $popsize --n-seq-indv 10 \
        --library-type ddRAD --enz PstI --enz2 EcoRI \
        --insert-mean 800 \
        --pcr-cycles 5 --coverage 20 --read-length 150 

    writeTruthDataFiles $OUTDIR

else
    echo "Directory already exists: $OUTDIR"
fi

# by default indels have a %1 probability

#
# Write dummy .fq.gz files from .fa.gz radinition output
#

# Convert .fa files to .fq with seqtk
RAD_READS_DIR=$(echo $OUTDIR | sed 's,$,/rad_reads,')
echo "writing fastq files with dummy quality scores (@)"
for fa in "${RAD_READS_DIR}"/*.fa.gz; do
    fq=$(echo $fa | sed 's/.fa.gz/.fq.gz/')
    if [ -f "$fq" ] && [ -s "$fq" ]; then
        echo "$fq exists and is not empty."
    else
        #echo "$fq does not exist or is empty."
        singularity exec seqtk%3A1.4--he4a0461_1 seqtk seq -F '@' $fa | gzip -c > $fq
    fi
done

#
# Write Samplesheet for nf-core/radseq
#
echo "Writing samplesheet: "${RAD_READS_DIR}/$DIRNAME.csv""
# Write SampleSheet
echo "sample,fastq_1,fastq_2,umi_barcodes" > "${RAD_READS_DIR}/$DIRNAME.csv"
paste -d',' <(for i in "${RAD_READS_DIR}"/*.1.fq.gz; do basename $i | cut -f1 -d'.' -; done)\
 <(ls "${RAD_READS_DIR}"/*.1.fq.gz) <(ls "${RAD_READS_DIR}"/*.2.fq.gz)\
 <(for i in "${RAD_READS_DIR}"/*.1.fq.gz; do echo "false"; done)\
 >> "${RAD_READS_DIR}/$DIRNAME.csv"
 
#
# Pass to nf-core/radseq
#

nextflow pull Gabriel-A-Barrett/radseq
nextflow pull Gabriel-A-Barrett/nf-vg-pipeline
######################################################################################################################
#                                              REFERENCE
######################################################################################################################
# Run RADseq and findMostAccurateVCF.py
#
NF_RADSEQ_REFERENCE_BED="${OUTDIR}/nf-radseq/reference/bwa-mem2/intervals/bedops_merge/*.bed" # FILE
RADSEQ_REF_VCF_DIR="${OUTDIR}/nf-radseq/reference/variant_calling" # DIRECTORY
RAD_REF_VCF_DIR="${RADSEQ_REF_VCF_DIR}/filter"

nextflow run Gabriel-A-Barrett/radseq -r dev --input "${RAD_READS_DIR}/$DIRNAME.csv" -c ~/config/unity.config --sequence_type 'PE'\
 --method 'reference' --genome $GENOME --outdir "${OUTDIR}/nf-radseq" -resume --save_reference_indices 'true'

# Write chromosome and position file across reference vcf's
echo "creating _chrompos.txt inside ${RAD_REF_VCF_DIR}"
for vcf in ${RAD_REF_VCF_DIR}/*.vcf.gz; do
    name=$(basename $vcf | cut -f 1,2 -d '.')
    singularity exec bcftools%3A1.17--haef29d1_0 bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $vcf > "${RAD_REF_VCF_DIR}/${name}_chrompos.txt"
done

BEST_REF_VCF=$(paste -d '' <(echo "${RAD_REF_VCF_DIR}/") <(python3 findMostAccurateVCF.py "${RAD_REF_VCF_DIR}" True))
echo "Constructing vg based on: $BEST_REF_VCF"

calculate_population_statistics $BEST_REF_VCF $TRUTH_NORM_VCF "nf-radseq-reference" $OUTDIR

#
# Define variables for VG
# 
REF_FAI="/project/pi_jpuritz_uri_edu/radinitio/${RADSEQ_REF_VCF_DIR}/../samtools/index/*.fasta.fai"
singularity exec htslib%3A1.18--h81da01d_0 tabix -p vcf -f $BEST_REF_VCF # radseq doesn't index reference
REF_TBI="${BEST_REF_VCF}.tbi"
FQ="${RAD_READS_DIR}/*.{1,2}.fq.gz"

echo "Constructing variant graph based on: $BEST_REF_VCF "

nextflow run Gabriel-A-Barrett/nf-vg-pipeline --fasta ${GENOME} --fai "${REF_FAI}"\
 --vcf "${BEST_REF_VCF}" --tbi "${REF_TBI}" -resume\
 --fq "${FQ}" --output_mode 'symlink' -c ~/config/unity.config --outdir "${OUTDIR}/vg/reference" --remove_duplicates true

RAD_VG_REF_VCF_DIR="${OUTDIR}/vg/reference/BCFTOOLS/NORM"
singularity exec bcftools%3A1.17--haef29d1_0 bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "${RAD_VG_REF_VCF_DIR}/msp_norm.vcf.gz" > "${RAD_VG_REF_VCF_DIR}/msp_norm_chrompos.txt"
BEST_VG_REF_VCF=$(paste -d '' <(echo "${RAD_VG_REF_VCF_DIR}/") <(python3 findMostAccurateVCF.py "${RAD_VG_REF_VCF_DIR}" True))
echo "${BEST_VG_REF_VCF}"
singularity exec htslib%3A1.18--h81da01d_0 tabix -f -p vcf "${BEST_VG_REF_VCF}"

# Function measures Fst and pi 
#   Runs on imputted genotypes
#   outputs: _{stat}.txt & _truth_{stat}.txt
calculate_population_statistics $BEST_VG_REF_VCF $TRUTH_NORM_VCF "nf-vg-reference" $OUTDIR

# Move files to ./test_accuracy
BEST_REF_VCF_CHROMPOS=$(echo "${BEST_REF_VCF}" | sed 's/.vcf.gz/_chrompos.txt/')
BEST_VG_REF_VCF_CHROMPOS=$(echo "${BEST_VG_REF_VCF}" | sed 's/.vcf.gz/_chrompos.txt/')
cp ${BEST_REF_VCF_CHROMPOS} "${OUTDIR}/test_accuracy/nf-radseq-reference_chrompos.txt"
cp ${BEST_VG_REF_VCF_CHROMPOS} "${OUTDIR}/test_accuracy/nf-vg-reference_chrompos.txt"

echo "Finished reference section"

######################################################################################################################
#                                              DENOVO 
######################################################################################################################
RADSEQ_DENOVO_VCF_DIR="${OUTDIR}/nf-radseq/denovo/variant_calling" # DIRECTORY
NF_RADSEQ_DENOVO_BED=$(echo $NF_RADSEQ_REFERENCE_BED | sed 's/reference/denovo/') # FILE
RAD_DENOVO_VCF_DIR=$(echo $RAD_REF_VCF_DIR | sed 's/reference/denovo/')

nextflow run Gabriel-A-Barrett/radseq -r dev --input "${RAD_READS_DIR}/$DIRNAME.csv" -c ~/config/unity.config --sequence_type 'PE'\
 --method 'denovo' --outdir "${OUTDIR}/nf-radseq" -resume --save_reference_indices 'true'\
 --minreaddepth_withinindividual=4 --minreaddepth_betweenindividual=5 --only_use_best_reference true
# Write chromosome and position file across reference vcf's
#module load bcftools/1.15
echo "creating _chrompos.txt inside $RAD_DENOVO_VCF_DIR"
for vcf in ${RAD_DENOVO_VCF_DIR}/*.vcf.gz; do
    name=$(basename $vcf | cut -f 1,2 -d '.')
    singularity exec bcftools%3A1.17--haef29d1_0 bcftools query -f '%CHROM\t%POS\n' $vcf > "${RAD_DENOVO_VCF_DIR}/${name}_chrompos.txt"
done
BEST_DENOVO_VCF=$(paste -d '' <(echo "${RAD_DENOVO_VCF_DIR}/") <(python3 findMostAccurateVCF.py "${RAD_DENOVO_VCF_DIR}/" True))
echo "best denovo vcf: $BEST_DENOVO_VCF"
calculate_population_statistics $BEST_DENOVO_VCF $TRUTH_NORM_VCF "nf-radseq-denovo" $OUTDIR
#
# Define variables for vg in denovo mode
#
DENOVO_FASTA=$(echo $BEST_DENOVO_VCF | sed -E 's/msp_([0-9]+)_([0-9]+)_.+\.vcf\.gz/msp_\1_\2_rainbow.fasta/' | sed 's,variant_calling/filter/,write_fasta/,') # take the denovo vcf outputted by BEST_DENOVO_VCF and use denovo combination to select reference
DENOVO_FAI=$(paste -d '' <(echo "$OUTDIR/nf-radseq/denovo/samtools/index/") <(basename $BEST_DENOVO_VCF | sed -E 's/msp_([0-9]+)_([0-9]+)_.+\.vcf\.gz/msp_\1_\2_rainbow.fasta/'))
singularity exec htslib%3A1.18--h81da01d_0 tabix -p vcf -f ${BEST_DENOVO_VCF} # radseq doesn't index reference
DENOVO_TBI="${BEST_DENOVO_VCF}.tbi"
FQ="${RAD_READS_DIR}/*.{1,2}.fq.gz"
echo "Variant Graph Mode: Denovo"
echo "   FASTA: $DENOVO_FASTA"
echo "   FAI: $DENOVO_FAI"
echo "   VCF: $BEST_DENOVO_VCF"
echo "   TBI: $DENOVO_TBI"
nextflow run Gabriel-A-Barrett/nf-vg-pipeline --fasta ${GENOME} --fai "${DENOVO_FAI}"\
 --vcf "${BEST_DENOVO_VCF}" --tbi "${DENOVO_TBI}"\
 --fq "${FQ}" -resume --output_mode 'symlink' -c ~/config/unity.config --outdir "${OUTDIR}/vg/denovo" --remove_duplicates true
RAD_VG_DENOVO_VCF_DIR=$(echo $RAD_VG_REF_VCF_DIR | sed 's/reference/denovo/')
singularity exec bcftools%3A1.17--haef29d1_0 bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "${RAD_VG_DENOVO_VCF_DIR}/msp_norm.vcf.gz" > "${RAD_VG_DENOVO_VCF_DIR}/msp_norm_chrompos.txt"
BEST_VG_DENOVO_VCF=$(paste -d '' <(echo "${RAD_VG_DENOVO_VCF_DIR}/") <(python3 findMostAccurateVCF.py "${RAD_VG_DENOVO_VCF_DIR}" True))
echo "${BEST_VG_DENOVO_VCF}"
singularity exec htslib%3A1.18--h81da01d_0 tabix -f -p vcf "${BEST_VG_DENOVO_VCF}"

# Function measures Fst and pi 
#   Runs on imputted genotypes
#   outputs: _{stat}.txt & _truth_{stat}.txt
calculate_population_statistics $BEST_VG_DENOVO_VCF $TRUTH_NORM_VCF "nf-vg-denovo" $OUTDIR



