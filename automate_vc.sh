#!/bin/bash

# stop the script as soon as an error occurs
set -euo pipefail

#----------------------------- Description -----------------------------
# This script performs variant calling using BCFtools.
# It generates a raw BCF file from alignment BAM files listed in a bam_list,
# calls variants with ploidy set to 1 (appropriate for haploid organisms),
# converts the BCF output to VCF format, compresses and indexes the final VCF,
# and generates variant statistics.
#
# Usage: bash automate_vc.sh
#
# Dependencies: bcftools, bgzip

#----------------------------dependency checks----------------------------------
for cmd in bcftools bgzip; do
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

# Define directories alignment
align_dir="alignment"
bam_list="$align_dir/bam_list.txt"
variants_dir="variants"

# Create necessary directories if they do not exist.
mkdir -p "$variants_dir" 

#--------------------------- Step 1: Generate raw BCF file ---------------------------
echo "========== Generating raw BCF file =========="
# Generate a BCF file from the list of BAM files using bcftools mpileup.
bcftools mpileup -Ou -f "$ref" -b "$bam_list" > "$variants_dir/variants.bcf"

#--------------------------- Step 2: Call Variants ---------------------------
echo "========== Calling variants with ploidy level set to 1 =========="
# Call variants from the raw BCF file with ploidy set to 1 (haploid).
bcftools call -mv --ploidy 1 -Ob -o "$variants_dir/called_variants.bcf" "$variants_dir/variants.bcf"

#--------------------------- Step 2: Convert BCF to VCF ------------------------
echo "========== Converting BCF to VCF format =========="
# Convert the BCF output to VCF format.
bcftools view "$variants_dir/called_variants.bcf" > "$variants_dir/called_variants.vcf"

#------------------- Step 3: Index final VCF and generate stats ----------------
echo "========== Compressing and indexing the final VCF file =========="
# Compress the VCF file with bgzip and then index it with bcftools.
bgzip -c "$variants_dir/called_variants.vcf" > "$variants_dir/called_variants.vcf.gz"
bcftools index "$variants_dir/called_variants.vcf.gz"

# Generate statistics on the final VCF file.
bcftools stats "$variants_dir/called_variants.vcf.gz" > "$variants_dir/bcf_stats.txt"

#--------------------------- Completion Message ---------------------------
echo "========== Variant calling completed successfully =========="
echo "Check VCF file and statistics in: $variants_dir"
