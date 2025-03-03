#!/bin/bash

# stop the script as soon as an error occurs
set -euo pipefail

#-------------------------------description-------------------------------------
# This script aligns trimmed reads to a reference genome using BWA.
# It indexes the reference genome, aligns each sample, converts SAM files to
# sorted BAM files, indexes the BAM files, and generates alignment statistics.
# It automatically detects samples in the "trimmed_reads" directory based on file
# naming conventions.
# Supported naming patterns are sample_R1.fastq.gz / sample_R2.fastq.gz or 
# sample_1.fastq.gz  / sample_2.fastq.gz
#
# Usage: bash automate_algnmnt.sh
#
# Dependencies: bwa, samtools

#---------------------------dependency checks---------------------------
# Verify that required commands (tools) are installed and available in systems PATH
for cmd in bwa samtools picard; do
    if ! command -v "$cmd" &> /dev/null; then
        echo "Error: $cmd is not installed or not in PATH."
        exit 1
    fi
done

#---------------------------directory setup-------------------------------------
# Automatically search for a FASTA (.fasta) or FNA (.fna) reference file in the "reference" directory.
ref=$(find reference -maxdepth 1 -type f \( -iname "*.fasta" -o -iname "*.fna" \) | head -n 1)
if [ -z "$ref" ]; then
    echo "Error: No reference file (.fasta or .fna) found in the 'reference' directory."
    exit 1
fi
echo "Using reference file: $ref"

# Define directories for trimmed reads, alignment outputs, and statistics.
trim_dir="trimmed_reads"
align_dir="alignment"
stats_dir="alignment/stats"

# Create necessary directories if they do not exist.
mkdir -p "$align_dir" "$stats_dir"

#---------------------------automate sample names detection---------------------
# Use an array to store detected sample names
samples=()

# Detect samples with the "_R1" pattern
for file in "$trim_dir"/*_R1*; do
    if [[ -f "$file" ]]; then
        base=$(basename "$file")
        # extract the sample name: everything before '_R1'
        sample="${base%%_R1*}"
        # add sampe name to samples array
        samples+=("$sample") 
    fi
done

# Detect samples with the "_1" pattern if samples array is still empty
if [ ${#samples[@]} -eq 0 ]; then
    for file in "$trim_dir"/*_1*; do
        if [[ -f "$file" ]]; then
            base=$(basename "$file")
            # Extract the sample name: everything before '_1'
            sample="${base%%_1*}"
            # add sampe name to samples array
            samples+=("$sample")
        fi
    done
fi

# Exit if no samples were detected.
if [ ${#samples[@]} -eq 0 ]; then
    echo "Error: No trimmed read files found in $trim_dir with expected naming pattern."
    exit 1
fi

#-----------------Step 1: Index the Reference Genome----------------------------
echo "==========Indexing reference genome: $ref=========="
bwa index "$ref"

#-----------------Steps 2&3: alignment and post-alignment processing------------
# Loop through each detected sample
for sample in "${samples[@]}"; do
    # Use the "_R1/_R2" naming pattern
    R1_files=( "$trim_dir/${sample}_R1"* )
    R2_files=( "$trim_dir/${sample}_R2"* )
    if [[ -f "${R1_files[0]}" && -f "${R2_files[0]}" ]]; then
        R1="${R1_files[0]}"
        R2="${R2_files[0]}"
    else
        # If not found, use the "_1/_2" naming pattern
        R1_files=( "$trim_dir/${sample}_1"* )
        R2_files=( "$trim_dir/${sample}_2"* )
        if [[ -f "${R1_files[0]}" && -f "${R2_files[0]}" ]]; then
            R1="${R1_files[0]}"
            R2="${R2_files[0]}"
        else
            echo "Warning: Paired files for sample '$sample' not found with expected patterns. Skipping."
            continue
        fi
    fi
    
    # Define file paths.
    sam_file="$align_dir/${sample}.sam"
    bam_file="$align_dir/${sample}.sorted.bam"
    dedup_bam="$align_dir/${sample}.sorted.dedup.bam"
    metrics_file="$stats_dir/${sample}.dedup.txt"
    
    # Align reads to the reference genome using bwa mem.
    echo "========= Aligning $sample to the reference genome ========="
    bwa mem -t 4 -R "@RG\tID:$sample\tSM:$sample\tPL:illumina" "$ref" "$R1" "$R2" > "$sam_file"
    
    # Convert SAM to coordinate-sorted BAM using samtools.
    echo "========= Converting SAM to sorted BAM for $sample ========="
    samtools view -S -b "$sam_file" | samtools sort -o "$bam_file"
    
    # Remove the intermediate SAM file.
    rm -f "$sam_file"
    
    # Mark duplicates using Picard.
    echo "========= Marking duplicates for $sample with Picard ========="
    picard MarkDuplicates -INPUT "$bam_file" -OUTPUT "$dedup_bam" -METRICS_FILE "$metrics_file"

    # Index the deduplicated BAM file.
    echo "========= Indexing deduplicated BAM file for $sample ========="
    samtools index "$dedup_bam"
    
    # Generate alignment statistics on the deduplicated BAM file.
    echo "========= Generating alignment statistics for $sample ========="
    samtools flagstat "$dedup_bam" > "$stats_dir/${sample}.flagstat"
    samtools stats "$dedup_bam" > "$stats_dir/${sample}.stats"
done

echo "=========Generating MultiQC report for alignment========="
# aggregate alignment statistics
multiqc -f -o $stats_dir -n multiqc_alignment_report.html $stats_dir

#---------------------------Create BAM List---------------------------
# Create a text file listing all deduplicated BAM files.
bam_list="$align_dir/bam_list.txt"
echo "========= Creating bam_list.txt ========="
ls "$align_dir"/*.sorted.dedup.bam > "$bam_list"
echo "bam_list.txt created at $bam_list"

#---------------------------Completion Message---------------------------
echo "========= Alignment completed ========="
echo "Check deduplicated BAM files in: $align_dir and alignment statistics in: $stats_dir"

    