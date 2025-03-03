#!/bin/bash
set -euo pipefail

#---------------------------------Description-----------------------------------
# End-to-End Variant Analysis Pipeline
#
# This script performs the following steps:
# 1. Quality control and trimming of raw sequencing reads.
# 2. Alignment of trimmed reads to a reference genome, duplicate marking,
#    and generation of alignment statistics.
# 3. Variant calling using bcftools.
# 4. Annotation and effect prediction using snpEff.
#
# Dependencies: fastqc, multiqc, trim_galore, bwa, samtools, picard, bcftools,
#               bgzip, snpEff, Java.
#
# Usage: bash end_to_end_variant_analysis.sh

#------------------------------Dependency Checks-----------------------------------
# Verify that required commands (tools) are installed and available in systems PATH
for cmd in fastqc multiqc trim_galore bwa samtools picard bcftools bgzip java snpEff; do
    if ! command -v "$cmd" &> /dev/null; then
        echo "Error: $cmd is not installed or not in PATH."
        exit 1
    fi
done

#--------------------------------Directory Setup------------------------------------
# Define directories for raw reads, trimmed reads, alignment, variants, and annotation.
raw_dir="raw_reads"
trim_dir="trimmed_reads"
qc_raw_dir="raw_reads/qc_raw"
qc_trim_dir="trimmed_reads/qc_trimmed"
align_dir="alignment"
stats_dir="$align_dir/stats"
variants_dir="variants"
annotation_dir="annotation"

# Create directories if they do not exist
mkdir -p "$raw_dir" "$qc_raw_dir" "$trim_dir" "$qc_trim_dir" "$align_dir" "$stats_dir" "$variants_dir" "$annotation_dir"


#-------------------------Step 1: Quality Control & Trimming-------------------------
echo "===================Performing quality checks on raw reads===================="
# Run FastQC on all regular files in raw_dir.
fastqc -o "$qc_raw_dir" $(find "$raw_dir" -maxdepth 1 -type f)

echo "===========Generating MultiQC Report for Raw Reads==========================="
multiqc -f -o "$qc_raw_dir" -n multiqc_raw_report.html "$qc_raw_dir"

# Detect sample names using "_R1" or "_1" patterns.
samples=()
for file in "$raw_dir"/*_R1*; do
    if [[ -f "$file" ]]; then
        base=$(basename "$file")
        sample="${base%%_R1*}"
        samples+=("$sample")
    fi
done
if [ ${#samples[@]} -eq 0 ]; then
    for file in "$raw_dir"/*_1*; do
        if [[ -f "$file" ]]; then
            base=$(basename "$file")
            sample="${base%%_1*}"
            samples+=("$sample")
        fi
    done
fi
if [ ${#samples[@]} -eq 0 ]; then
    echo "Error: No raw read files found in $raw_dir with expected naming patterns (_R1/_R2 or _1/_2)."
    exit 1
fi


echo "==========Trimming reads and performing quality checks on trimmed reads=========="
for sample in "${samples[@]}"; do
    # Try _R1/_R2 pattern; if not, try _1/_2.
    R1_files=( "$raw_dir/${sample}_R1"* )
    R2_files=( "$raw_dir/${sample}_R2"* )
    if [[ -f "${R1_files[0]}" && -f "${R2_files[0]}" ]]; then
        R1="${R1_files[0]}"
        R2="${R2_files[0]}"
    else
        R1_files=( "$raw_dir/${sample}_1"* )
        R2_files=( "$raw_dir/${sample}_2"* )
        if [[ -f "${R1_files[0]}" && -f "${R2_files[0]}" ]]; then
            R1="${R1_files[0]}"
            R2="${R2_files[0]}"
        else
            echo "Warning: Paired files for sample '$sample' not found. Skipping."
            continue
        fi
    fi

    echo "Trimming sample: $sample"
    trim_galore --paired \
                --quality 20 \
                --fastqc \
                --fastqc_args "-o ${qc_trim_dir}" \
                --output_dir "${trim_dir}" \
                "$R1" "$R2"
done

echo "========== Generating MultiQC Reports for Trimmed Reads =========="
mv "$trim_dir"/*_trimming_report.txt "$qc_trim_dir"
multiqc -f -o "$qc_trim_dir" -n trim_galore_report.html "$qc_trim_dir"
multiqc -f -o "$qc_trim_dir" -n multiqc_trimmed_report.html "$qc_trim_dir"

echo "Quality control and trimming completed."
echo "Inspect raw reads report: $qc_raw_dir/multiqc_raw_report.html"
echo "Inspect trimmed reads report: $qc_trim_dir/multiqc_trimmed_report.html"
echo "Inspect trim-galore report: $qc_trim_dir/trim_galore_report.html"

#-------------------Step 2: Alignment & Post-Alignment Processing-----------------
# Automatically detect a reference file in FASTA (.fasta or .fna) format.
ref=$(find reference -maxdepth 1 -type f \( -iname "*.fasta" -o -iname "*.fna" \) | head -n 1)
if [ -z "$ref" ]; then
    echo "Error: No reference file (.fasta or .fna) found in the 'reference' directory."
    exit 1
fi
echo "Using reference file: $ref"

echo "=========================Indexing Reference Genome========================"
bwa index "$ref"

# trimmed reads
for sample in "${samples[@]}"; do
    # Use _R1/_R2 or _1/_2 pattern for trimmed reads.
    R1_files=( "$trim_dir/${sample}_1_val_1."* )
    R2_files=( "$trim_dir/${sample}_2_val_2."* )
    if [[ -f "${R1_files[0]}" && -f "${R2_files[0]}" ]]; then
        R1="${R1_files[0]}"
        R2="${R2_files[0]}"
    else
        R1_files=( "$trim_dir/${sample}_R1"* )
        R2_files=( "$trim_dir/${sample}_R2"* )
        if [[ -f "${R1_files[0]}" && -f "${R2_files[0]}" ]]; then
            R1="${R1_files[0]}"
            R2="${R2_files[0]}"
        else
            echo "Warning: Paired trimmed files for sample '$sample' not found. Skipping."
            continue
        fi
    fi

    sam_file="$align_dir/${sample}.sam"
    bam_file="$align_dir/${sample}.sorted.bam"
    dedup_bam="$align_dir/${sample}.sorted.dedup.bam"
    metrics_file="$stats_dir/${sample}.dedup.txt"

    echo "=================Aligning sample $sample==========================="
    bwa mem -t 4 -R "@RG\tID:$sample\tSM:$sample\tPL:illumina" "$ref" "$R1" "$R2" > "$sam_file"

    echo "========== Converting SAM to Sorted BAM for $sample =========="
    samtools view -S -b "$sam_file" | samtools sort -o "$bam_file"
    rm -f "$sam_file"

    echo "========== Marking Duplicates for $sample using Picard =========="
    picard MarkDuplicates -INPUT "$bam_file" -OUTPUT "$dedup_bam" -METRICS_FILE "$metrics_file"

    echo "========== Indexing Deduplicated BAM for $sample =========="
    samtools index "$dedup_bam"

    echo "========== Generating Alignment Statistics for $sample =========="
    samtools flagstat "$dedup_bam" > "$stats_dir/${sample}.flagstat"
    samtools stats "$dedup_bam" > "$stats_dir/${sample}.stats"
done

echo "========== Generating MultiQC Report for Alignment =========="
multiqc -f -o "$stats_dir" -n multiqc_alignment_report.html "$stats_dir"

# Create a list of deduplicated BAM files.
bam_list="$align_dir/bam_list.txt"
echo "========== Creating BAM List =========="
ls "$align_dir"/*.sorted.dedup.bam > "$bam_list"
echo "BAM list created at: $bam_list"

#-----------------------Step 3: Variant Calling---------------------------------
echo "==================================Variant Calling========================="

# Use the same reference FASTA file as for alignment
echo "Using reference file: $ref"

echo "====================Generating raw BCF file==============================="
bcftools mpileup -Ou -f "$ref" -b "$bam_list" > "$variants_dir/variants.bcf"

echo "========== Calling Variants with ploidy set to 1=========================="
bcftools call -mv --ploidy 1 -Ob -o "$variants_dir/called_variants.bcf" "$variants_dir/variants.bcf"

echo "========================Converting BCF to VCF============================="
bcftools view "$variants_dir/called_variants.bcf" > "$variants_dir/called_variants.vcf"

echo "========================Compressing and Indexing VCF======================"
bgzip -c "$variants_dir/called_variants.vcf" > "$variants_dir/called_variants.vcf.gz"
bcftools index "$variants_dir/called_variants.vcf.gz"

echo "===========================Generating Variant Statistics=================="
bcftools stats "$variants_dir/called_variants.vcf.gz" > "$variants_dir/bcf_stats.txt"

echo "Variant calling completed. Check VCF and stats in: $variants_dir"

#-----------------Step 4: Annotation and Effect Prediction with snpEff---------
echo "=========================Annotation with snpEff==========================="
# Dynamically detect a GenBank (.gb) reference file in the "reference" directory.
ref_gb=$(find reference -maxdepth 1 -type f -name "*.gb" | head -n 1)
if [ -z "$ref_gb" ]; then
    echo "Error: No GenBank reference file (.gb) found in the 'reference' directory."
    exit 1
fi
echo "Using reference GenBank file: $ref_gb"

# Extract organism details from the GenBank filename.
ref_base=$(basename "$ref_gb" .gb)
organ_db="$ref_base"
organism="${ref_base//_/ }"
genome="${organism}.genome"

# Dynamically detect snpEff configuration file from the Conda packages directory.
snpeff_config=$(find "$(dirname "$(dirname "$CONDA_PREFIX")")/pkgs" -type f -name "snpEff.config" | head -n 1)
if [ -z "$snpeff_config" ]; then
    echo "Error: Could not locate snpEff configuration file."
    exit 1
fi
echo "Using snpEff configuration file: $snpeff_config"

# Create directory structure for snpEff annotation
mkdir -p "$annotation_dir/data/$organ_db"

echo "======================Building Custom snpEff Database====================="
# Copy the GenBank file into the snpEff data directory
cp "$ref_gb" "$annotation_dir/data/$organ_db/genes.gbk"

# Copy the snpEff configuration file for local modification
cp "$snpeff_config" "$annotation_dir/snpeff_copy.config"

# Append organism and genome information to the copied configuration.
echo "# $organism" >> "$annotation_dir/snpeff_copy.config"
echo "$genome : $organism, complete genome" >> "$annotation_dir/snpeff_copy.config"

# Build the custom snpEff database.
snpEff build -genbank -c "$annotation_dir/snpeff_copy.config" -v "$organ_db"

echo "=============================Annotating Variants with snpEf================"
# Annotate variants; output annotated VCF and stats.
snpEff ann -o vcf \
           -stats "$annotation_dir/snpEff_summary.html" \
           -c "$annotation_dir/snpeff_copy.config" \
           "$organ_db" \
           "$variants_dir/called_variants.vcf.gz" > "$annotation_dir/ann_called_variants.vcf"

echo "Annotation and effect prediction completed."
echo "Check annotated VCF in $annotation_dir/ann_called_variants.vcf"
echo "Check snpEff summary files: in $annotation_dir"
echo "=============Variant Analysis Completed==============="
