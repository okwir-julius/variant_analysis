#!/bin/bash

# stop the script as soon as an error occurs
set -euo pipefail

#------------------------------description--------------------------------------
# This script processes sequencing data by performing quality checks, trimming 
# low-quality reads, and re-checking quality post-trimming.
# It automatically detects samples from the "raw_reads" directory based on file
# naming conventions. 
# Supported naming patterns are _R1/_R2 or _1/_2
#
# Usage: bash automate_qc.sh
#
# Dependencies: fastqc, multiqc, trim_galore

#---------------------------dependency checks-----------------------------------
# Verify that required commands (tools) are installed and available in systems PATH
for cmd in fastqc multiqc trim_galore; do
    if ! command -v "$cmd" &> /dev/null; then
        echo "Error: $cmd is not installed or not in PATH."
        exit 1
    fi
done

#---------------------------directory setup-------------------------------------
# Define directories for raw reads, trimmed reads, and QC outputs
raw_dir="raw_reads"
trim_dir="trimmed_reads"
qc_raw_dir="raw_reads/qc_raw"
qc_trim_dir="trimmed_reads/qc_trimmed"

# Create directories if absent.
mkdir -p "$raw_dir" "$qc_raw_dir" "$trim_dir" "$qc_trim_dir"

#---------------------------automate sample names detection---------------------
# Use an array to store detected sample names
samples=()

# Detect samples with the "_R1" pattern
for file in "$raw_dir"/*_R1*; do
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
    for file in "$raw_dir"/*_1*; do
        if [[ -f "$file" ]]; then
            base=$(basename "$file")
            # Extract the sample name: everything before '_1'
            sample="${base%%_1*}"
            # add sampe name to samples array
            samples+=("$sample")
        fi
    done
fi

# If still no samples are found, exit with an error
if [ ${#samples[@]} -eq 0 ]; then
    echo "Error: No raw read files found in $raw_dir with expected naming patterns, \n
    e.g _R1/_R2 or _1/_2"
    exit 1
fi

#---------------------------Step 1: quality check on raw reads------------------
echo "=========Performing quality checks on raw reads========="
# Use 'find' to select only regular files in the raw_dir (ignore subdirectories)
# and run FastQC on them
fastqc -o "$qc_raw_dir" $(find "$raw_dir" -maxdepth 1 -type f)

echo "=========Generating MultiQC report for raw reads========="
# Generate a MultiQC report from the FastQC output. -f overwrites existing report
multiqc -f -o "$qc_raw_dir" -n multiqc_raw_report.html "$qc_raw_dir"

#------------------Steps 2&3: trim reads and check quality of trimmed reads-----
echo "======Trimming reads and performing quality checks on trimmed reads======="
# Loop through each detected sample
for sample in "${samples[@]}"; do
    # Use the "_R1/_R2" naming pattern
    R1_files=( "$raw_dir/${sample}_R1"* )
    R2_files=( "$raw_dir/${sample}_R2"* )
    if [[ -f "${R1_files[0]}" && -f "${R2_files[0]}" ]]; then
        R1="${R1_files[0]}"
        R2="${R2_files[0]}"
    else
        # If not found, use the "_1/_2" naming pattern
        R1_files=( "$raw_dir/${sample}_1"* )
        R2_files=( "$raw_dir/${sample}_2"* )
        if [[ -f "${R1_files[0]}" && -f "${R2_files[0]}" ]]; then
            R1="${R1_files[0]}"
            R2="${R2_files[0]}"
        else
            echo "Warning: Paired files for sample '$sample' not found with expected patterns. Skipping."
            continue
        fi
    fi

    echo "Trimming and running FastQC on sample: $sample"
    # trim, run FastQC on trimmed output, and redirect fastqc output to $qc_trim_dir
    trim_galore --paired \
                --quality 20 \
                --fastqc \
                --fastqc_args "-o ${qc_trim_dir}" \
                --output_dir "${trim_dir}" \
                "${R1}" "${R2}"
    
done

echo "=========Generating MultiQC report for trim-galore========="
# move trim galore report from $trim_dir to $qc_trim_dir
mv "$trim_dir"/*_trimming_report.txt "$qc_trim_dir"

# generate report
multiqc -f -o "$qc_trim_dir" -n trim_galore_report.html "$qc_trim_dir"

echo "=========Generating MultiQC report for trimmed reads========="
# Generate a MultiQC report for the FastQC output of the trimmed reads
multiqc -f -o "$qc_trim_dir" -n multiqc_trimmed_report.html "$qc_trim_dir"

#---------------------------Completion message---------------------------
echo "=========Quality control completed========="
echo "Inspect raw reads report in: $qc_raw_dir/multiqc_raw_report.html."
echo "Inspect trimmed reads report in: $qc_trim_dir/multiqc_trimmed_report.html."
echo "Inspect trim-galore report in: $qc_trim_dir/trim_galore_report.html."
