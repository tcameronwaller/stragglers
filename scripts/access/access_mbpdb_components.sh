#!/bin/bash

#chmod u+x script.sh

################################################################################
# Notes:
# first draft of script: TCW; 31 August 2022
# last revision of script: TCW; 1 September 2022

################################################################################
# Organize paths to files and directories.

# Read private, local file path for main project directory.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_storage=$(<"./storage_waller_metabolism.txt")
path_directory_mayo_bipolar_phenotype=$(<"./mayo_bipolar_phenotypes.txt")
path_directory_mayo_bipolar_genotype=$(<"./mayo_bipolar_genotype.txt")

# Define paths.
name_file_identifier="210421_id_matching_gwas.csv" # Identifier matching for phenotype records to GWAS1 and GWAS2 genotype records
name_file_phenotype="220513_BP_phenotypes.csv" # Richard S. Pendegraft prepared and shared this file on 13 May 2022
name_file_genetic_sex_case="MERGED.maf0.01.dosR20.8.noDups.fam" # Reference for genetic sex and case status
path_file_identifier_source="${path_directory_mayo_bipolar_phenotype}/${name_file_identifier}"
path_file_phenotype_source="${path_directory_mayo_bipolar_phenotype}/${name_file_phenotype}"
path_file_genetic_sex_case_source="${path_directory_mayo_bipolar_genotype}/${name_file_genetic_sex_case}"

path_directory_persistence="${path_directory_storage}/phenotypes/mayo_clinic_bipolar_disorder_access_2022-09-01"
path_file_identifier_persistence="${path_directory_persistence}/${name_file_identifier}"
path_file_phenotype_persistence="${path_directory_persistence}/${name_file_phenotype}"
path_file_genetic_sex_case_persistence="${path_directory_persistence}/${name_file_genetic_sex_case}"

path_directory_dock="${path_directory_process}/dock"
path_directory_access="${path_directory_dock}/access/mayo_bpdb"
path_file_identifier_product="${path_directory_access}/${name_file_identifier}"
path_file_phenotype_product="${path_directory_access}/${name_file_phenotype}"
path_file_genetic_sex_case_product="${path_directory_access}/${name_file_genetic_sex_case}"

# Echo each command to console.
set -x

# Remove previous versions.
rm -r $path_directory_access

# Create directory.
mkdir -p $path_directory_access

################################################################################
# Copy files to persistence directory.

if true; then
  # Remove previous versions.
  rm -r $path_directory_persistence
  # Create directory.
  mkdir -p $path_directory_persistence
  # Copy.
  cp $path_file_identifier_source $path_file_identifier_persistence
  cp $path_file_phenotype_source $path_file_phenotype_persistence
  cp $path_file_genetic_sex_case_source $path_file_genetic_sex_case_persistence
fi

################################################################################
# Copy files to processing directory.

cp $path_file_identifier_source $path_file_identifier_product
cp $path_file_phenotype_source $path_file_phenotype_product
cp $path_file_genetic_sex_case_source $path_file_genetic_sex_case_product
