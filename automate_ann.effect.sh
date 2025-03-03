#!/bin/bash

# stop the script as soon as an error occurs
set -euo pipefail

#----------------------------- Description -----------------------------
# This script builds a custom snpEff database and annotates variants.
# It uses a reference GenBank file to automatically derive organism details,
# constructs a custom snpEff database, modifies the snpEff configuration file,
# and then annotates variants in a VCF file to predict their effects.
#
# Usage: bash annotate_variants.sh
#
# Dependencies: snpEff, Java

#---------------------------dependency checks ---------------------------
# Check for Java, which is required to run snpEff.
for cmd in java snpEff; do
    if ! command -v "$cmd" &> /dev/null; then
        echo "Error: $cmd is not installed or not in PATH."
        exit 1
    fi
done

#---------------------------Files and directory setup ---------------------------
# Automatically search for a GenBank (.gb) reference file in the "reference" directory.
ref_gb=$(find reference -maxdepth 1 -type f -name "*.gb" | head -n 1)
if [ -z "$ref_gb" ]; then
    echo "Error: No GenBank reference file (.gb) found in the 'reference' directory."
    exit 1
fi
echo "Using reference GenBank file: $ref_gb"

# Automatically extract organism details from the ref_gb filename.
# For example, if ref_gb is "Salmonella_enterica.gb":
#   - ref_base becomes "Salmonella_enterica"
#   - organ_db is set to "Salmonella_enterica"
#   - organism is converted to "Salmonella enterica"
#   - genome is set to "Salmonella enterica.genome"
ref_base=$(basename "$ref_gb" .gb)
organ_db="$ref_base"
organism="${ref_base//_/ }"
genome="${organism}.genome"

# automatically detect the snpEff configuration file using the Conda packages directory
# Search the conda pkgs directory relative to CONDA_PREFIX
snpeff_config=$(find "$(dirname "$(dirname "$CONDA_PREFIX")")/pkgs" -type f -name "snpEff.config" | head -n 1)
if [ -z "$snpeff_config" ]; then
    echo "Error: Could not locate snpEff configuration file."
    exit 1
fi
echo "Using snpEff configuration file: $snpeff_config"

# set path to config file manually in case of error
# snpeff_config="path/to/snpEff.config" 

# Define other directories and configuration.
variants_dir="variants"                   
annotation_dir="annotation"                   

# Create necessary output directories for annotation.
mkdir -p "$annotation_dir/data/$organ_db"

#--------------------------- Build Custom snpEff Database ---------------------------
# Copy the GenBank file into the snpEff data directory for the organism.
cp "$ref_gb" "$annotation_dir/data/$organ_db/genes.gbk"

# Copy the snpEff configuration file to the annotation directory for modification.
cp "$snpeff_config" "$annotation_dir/snpeff_copy.config"

# Append organism and genome information to the copied snpEff configuration file.
echo "# $organism" >> "$annotation_dir/snpeff_copy.config"
echo "$genome : $organism, complete genome" >> "$annotation_dir/snpeff_copy.config"

# Build the custom snpEff database using the GenBank file.
snpEff build -genbank -c "$annotation_dir/snpeff_copy.config" -v "$organ_db"

#--------------------------- Annotation and Effect Prediction ---------------------------
# Annotate variants and predict their effects using the custom database.
# The final VCF is written to "ann_called_variants.vcf" in the annotation directory
snpEff ann -o vcf \
           -stats "$annotation_dir/snpEff_summary.html" \
           -c "$annotation_dir/snpeff_copy.config" "$organ_db" "$variants_dir/called_variants.vcf.gz" > \
           "$annotation_dir/ann_called_variants.vcf"


#--------------------------- Completion Message ---------------------------
echo "Annotation and effect prediction completed successfully."
echo "Check final VCF file in $annotation_dir/ann_called_variants.vcf"
echo "Check snpEff and Genes summary files in $annotation_dir directory"





