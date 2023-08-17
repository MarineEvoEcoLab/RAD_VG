#!/usr/bin/env bash
module load bcftools/1.15
module load htslib/1.12

# Check number of arguments is atleast one

if [ $# -ne 1 ]; then
    echo "Incorrect number of arguments. Please provide exactly one argument."
    exit 1
fi
echo "Argument provided: $1"
GENOME=$1
CHROMOSOME_LIST=$(echo $GENOME | sed 's/\(\.fa\|\.fasta\)$/_chromlist.txt/')
# Default outdir
OUTDIR=$(echo $GENOME | sed 's/\(\.fa\|\.fasta\)$//') # DIRECTORY
DIRNAME=$(basename $OUTDIR) # PREFIX

TRUTH_VCF="${OUTDIR}/ref_loci_vars/ri_master.vcf.gz" # FILE
TRUTH_NORM_VCF=$(echo $TRUTH_VCF | sed 's/ri_master/ri_master_norm/') # FILE

NF_RADSEQ_REFERENCE_BED=$OUTDIR/nf-radseq/reference/bwa-mem2/intervals/bedops_merge/msp.bed # FILE
NF_RADSEQ_DENOVO_BED=$(echo $NF_RADSEQ_REFERENCE_BED | sed 's/reference/denovo/') # FILE
RADSEQ_REF_VCF_DIR="${OUTDIR}/nf-radseq/reference/variant_calling" # DIRECTORY

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
        module -q load seqtk/2023-05-11
        seqtk seq -F '@' $fa | gzip -c > $fq
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
module load nextflow/23.04.1
module load singularity/3.7.0

nextflow pull Gabriel-A-Barrett/radseq
nextflow pull Gabriel-A-Barrett/nf-vg-pipeline
##############
# REFERENCE
##############
# radseq
nextflow run Gabriel-A-Barrett/radseq -r dev --input "${RAD_READS_DIR}/$DIRNAME.csv" -c ~/config/unity.config --sequence_type 'PE'\
 --method 'reference' --genome $GENOME --outdir "${OUTDIR}/nf-radseq" -resume --save_reference_indices true

mkdir -p $OUTDIR/test_accuracy
bcftools query -R $NF_RADSEQ_REFERENCE_BED -f '%CHROM\t%POS\n' $TRUTH_NORM_VCF > $OUTDIR/test_accuracy/ri_master_norm_nf_ref_chrom_pos.txt
# RAW VCF
bcftools query -f '%CHROM\t%POS\n' ${RADSEQ_REF_VCF_DIR}/msp_norm.vcf.gz > $OUTDIR/test_accuracy/msp_norm.vcf.chrom.pos.txt
# SMALL FILTER
bcftools query -f '%CHROM\t%POS\n' ${RADSEQ_REF_VCF_DIR}/filter/msp.vcf.gz > $OUTDIR/test_accuracy/msp_norm_filter.vcf.chrom.pos.txt
# HARD VCF
bcftools query -f '%CHROM\t%POS\n' ${RADSEQ_REF_VCF_DIR}/filter/msp_RADseqPopFilters.vcf.gz > $OUTDIR/test_accuracy/msp_RADseqPopFilters.vcf.chrom.pos.txt
# FILTERED VCF W/ rmved indv
#bcftools query -f '%CHROM\t%POS\n' ${RADSEQ_REF_VCF_DIR}/filter/msp_rmvindv_RADseqPopFilters.bcf.gz > $OUTDIR/test_accuracy/msp_rmvindv_RADseqPopFilters.vcf.chrom.pos.txt

python3 measureVcfAccuracy.py $OUTDIR/test_accuracy/ri_master_norm_nf_ref_chrom_pos.txt $OUTDIR/test_accuracy/msp_norm.vcf.chrom.pos.txt > $OUTDIR/test_accuracy/msp_norm.vcf.accuracy.txt
python3 measureVcfAccuracy.py $OUTDIR/test_accuracy/ri_master_norm_nf_ref_chrom_pos.txt $OUTDIR/test_accuracy/msp_norm_filter.vcf.chrom.pos.txt > $OUTDIR/test_accuracy/msp_norm_filter.vcf.accuracy.txt
python3 measureVcfAccuracy.py $OUTDIR/test_accuracy/ri_master_norm_nf_ref_chrom_pos.txt $OUTDIR/test_accuracy/msp_RADseqPopFilters.vcf.chrom.pos.txt > $OUTDIR/test_accuracy/msp_RADseqPopFilters.vcf.accuracy.txt
#python3 measureVcfAccuracy.py $OUTDIR/test_accuracy/ri_master_norm_nf_ref_chrom_pos.txt $OUTDIR/test_accuracy/msp_rmvindv_RADseqPopFilters.vcf.chrom.pos.txt

#echo "creating variant graph and calling genotypes"
# vg

REF_FAI="${RADSEQ_REF_VCF_DIR}/../samtools/index/${DIRNAME}.fasta.fai"
REF_VCF="${RADSEQ_REF_VCF_DIR}/filter/msp.vcf.gz"
tabix -p vcf $REF_VCF # radseq doesn't index reference
REF_TBI="/project/pi_jpuritz_uri_edu/radinitio/${RADSEQ_REF_VCF_DIR}/filter/msp.vcf.gz.tbi"
FQ="${RAD_READS_DIR}/*.{1,2}.fq.gz"

nextflow run Gabriel-A-Barrett/nf-vg-pipeline --fasta ${GENOME} --fai "${REF_FAI}"\
 --vcf "${REF_VCF}" --tbi "${REF_TBI}"\
 --fq "${FQ}" -resume --output_mode 'symlink' -c ~/config/unity.config --outdir "${OUTDIR}/vg"

bcftools query -f '%CHROM\t%POS\n' $OUTDIR/vg/BCFTOOLS/msp.vcf.gz > $OUTDIR/test_accuracy/msp_vg.vcf.chrom.pos.txt
python3 measureVcfAccuracy.py $OUTDIR/test_accuracy/ri_master_norm_nf_ref_chrom_pos.txt $OUTDIR/test_accuracy/msp_vg.vcf.chrom.pos.txt > $OUTDIR/test_accuracy/msp_vg.vcf.accuracy.txt

##########
# DENOVO
##########

nextflow run Gabriel-A-Barrett/radseq -r dev --input "${RAD_READS_DIR}/$DIRNAME.csv" -c ~/config/unity.config --sequence_type 'PE'\
 --method 'denovo' --outdir "${OUTDIR}/nf-radseq" -resume --save_reference_indices 'true'

echo "finished radseq denovo"

#
# Pass nf-core/radseq to vg
#

#module load bcftools/1.15

# Re-annotate chromosome Field in VCF
#paste <(bcftools query -f '%CHROM\n' "${OUTDIR}/msprime_vcfs/Chr1.vcf.bgz" | sort -n | uniq) <(cat $CHROMOSOME_LIST) > "${OUTDIR}/msprime_vcfs/reannotate.txt"

#bcftools annotate --rename-chrs "${OUTDIR}/msprime_vcfs/reannotate.txt" -O b -o "${OUTDIR}/msprime_vcfs/Chr1_rename.vcf.bgz" "${OUTDIR}/msprime_vcfs/Chr1.vcf.bgz"

#singularity exec vg.sif vg construct -r $GENOME -v $NF_RADSEQ_VCF > graph.vg

#singularity exec vg.sif vg index -x graph.xg -g graph.gcsa graph.vg

#singularity exec vg.sif vg autoindex --workflow giraffe --prefix /path/to --ref-fasta reference.fasta --vcf variants.vcf.gz

#singularity exec vg.sif vg giraffe --threads 10 -Z girraffe.gbz -m pangenome.min -d pangenome.dist -f sim.1.fq -f sim.2.fq > mapped.gam

#singularity exec vg.sif vg stats -a mapped.gam


# multiple steps for variant calling

# compute read support, mapq >=5
#vg pack -x graph.vg -g mapped.gam -o aln.pack -Q 5
# vg snarls graph.vg > graph.snarls # recommended for larger files
#vg call x.xg -k aln.pack -s <sample> > 1.vcf 

# merge outputs from vg call
#bcftools merge -m all


#GENOME_CHR=$(cat $CHROMOSOME_LIST)
#MSPRIME_VCF_CHR=$(bcftools query -f '%CHROM\n' "${OUTDIR}/msprime_vcfs/Chr1.vcf.bgz" | sort -n | uniq)

