#!/bin/bash

#chmod u+x script.sh

################################################################################
# Notes:
# first draft of script: TCW; 31 August 2022
# last revision of script: TCW; 13 September 2022

################################################################################
# Organize paths to files and directories.

# Read private, local file path for main project directory.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_storage=$(<"./storage_waller_metabolism.txt")
path_directory_mayo_bipolar_phenotype=$(<"./mayo_bipolar_phenotypes.txt")
path_directory_mayo_bipolar_genotype=$(<"./mayo_bipolar_genotype.txt")
path_directory_mayo_bipolar_genotype_pca=$(<"./mayo_bipolar_genotype_pca.txt")

# Define paths.

name_file_identifier="210421_id_matching_gwas.csv" # file date: 21 April 2021
name_file_phenotype="220513_BP_phenotypes.csv" # file date: 13 May 2022; Richard S. Pendegraft prepared and shared this file on 13 May 2022
name_file_genetic_sex_case="MERGED.maf0.01.dosR20.8.noDups.fam" # file date: 20 May 2022
name_file_genotype_pca="Top20_PCs.csv" # file date: 23 May 2022

path_file_identifier_source="${path_directory_mayo_bipolar_phenotype}/${name_file_identifier}"
path_file_phenotype_source="${path_directory_mayo_bipolar_phenotype}/${name_file_phenotype}"
path_file_genetic_sex_case_source="${path_directory_mayo_bipolar_genotype}/${name_file_genetic_sex_case}"
path_file_genotype_pca_source="${path_directory_mayo_bipolar_genotype_pca}/${name_file_genotype_pca}"

path_directory_persistence="${path_directory_storage}/phenotypes/mayo_clinic_bipolar_disorder_access_2022-09-13"
path_file_identifier_persistence="${path_directory_persistence}/${name_file_identifier}"
path_file_phenotype_persistence="${path_directory_persistence}/${name_file_phenotype}"
path_file_genetic_sex_case_persistence="${path_directory_persistence}/${name_file_genetic_sex_case}"
path_file_genotype_pca_persistence="${path_directory_persistence}/${name_file_genotype_pca}"

path_directory_dock="${path_directory_process}/dock"
path_directory_access="${path_directory_dock}/access/mayo_bpdb"
path_file_identifier_product="${path_directory_access}/${name_file_identifier}"
path_file_phenotype_product="${path_directory_access}/${name_file_phenotype}"
path_file_genetic_sex_case_product="${path_directory_access}/${name_file_genetic_sex_case}"
path_file_genotype_pca_product="${path_directory_access}/${name_file_genotype_pca}"

name_file_supplement_1="220325_BP_phenotypes.csv" # file date: 31 March 2021
name_file_supplement_2="211221_BP_phenotypes.csv" # file date: 23 December 2021
name_file_supplement_3="210902_BP_phenotypes.csv" # file date: 2 September 2021
name_file_supplement_4="210609_BP_phenotypes.csv" # file date: 9 June 2021
name_file_supplement_5="210422_BP_phenotypes.csv" # file date: 5 May 2021
name_file_supplement_6="210330_BP_phenotypes.csv" # file date: 13 April 2021
name_file_supplement_7="phen_bp.csv" # file date: 5 September 2022
path_file_supplement_1_source="${path_directory_mayo_bipolar_phenotype}/${name_file_supplement_1}"
path_file_supplement_2_source="${path_directory_mayo_bipolar_phenotype}/${name_file_supplement_2}"
path_file_supplement_3_source="${path_directory_mayo_bipolar_phenotype}/${name_file_supplement_3}"
path_file_supplement_4_source="${path_directory_mayo_bipolar_phenotype}/${name_file_supplement_4}"
path_file_supplement_5_source="${path_directory_mayo_bipolar_phenotype}/${name_file_supplement_5}"
path_file_supplement_6_source="${path_directory_mayo_bipolar_phenotype}/${name_file_supplement_6}"
path_file_supplement_7_source="${path_directory_mayo_bipolar_phenotype}/${name_file_supplement_7}"

path_file_supplement_1_persistence="${path_directory_persistence}/${name_file_supplement_1}"
path_file_supplement_2_persistence="${path_directory_persistence}/${name_file_supplement_2}"
path_file_supplement_3_persistence="${path_directory_persistence}/${name_file_supplement_3}"
path_file_supplement_4_persistence="${path_directory_persistence}/${name_file_supplement_4}"
path_file_supplement_5_persistence="${path_directory_persistence}/${name_file_supplement_5}"
path_file_supplement_6_persistence="${path_directory_persistence}/${name_file_supplement_6}"
path_file_supplement_7_persistence="${path_directory_persistence}/${name_file_supplement_7}"

path_file_supplement_1_product="${path_directory_access}/${name_file_supplement_1}"
path_file_supplement_2_product="${path_directory_access}/${name_file_supplement_2}"
path_file_supplement_3_product="${path_directory_access}/${name_file_supplement_3}"
path_file_supplement_4_product="${path_directory_access}/${name_file_supplement_4}"
path_file_supplement_5_product="${path_directory_access}/${name_file_supplement_5}"
path_file_supplement_6_product="${path_directory_access}/${name_file_supplement_6}"
path_file_supplement_7_product="${path_directory_access}/${name_file_supplement_7}"

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
  cp $path_file_genotype_pca_source $path_file_genotype_pca_persistence

  cp $path_file_supplement_1_source $path_file_supplement_1_persistence
  cp $path_file_supplement_2_source $path_file_supplement_2_persistence
  cp $path_file_supplement_3_source $path_file_supplement_3_persistence
  cp $path_file_supplement_4_source $path_file_supplement_4_persistence
  cp $path_file_supplement_5_source $path_file_supplement_5_persistence
  cp $path_file_supplement_6_source $path_file_supplement_6_persistence
  cp $path_file_supplement_7_source $path_file_supplement_7_persistence
fi

################################################################################
# Copy files to processing directory.

cp $path_file_identifier_source $path_file_identifier_product
cp $path_file_phenotype_source $path_file_phenotype_product
cp $path_file_genetic_sex_case_source $path_file_genetic_sex_case_product
cp $path_file_genotype_pca_source $path_file_genotype_pca_product

cp $path_file_supplement_1_source $path_file_supplement_1_product
cp $path_file_supplement_2_source $path_file_supplement_2_product
cp $path_file_supplement_3_source $path_file_supplement_3_product
cp $path_file_supplement_4_source $path_file_supplement_4_product
cp $path_file_supplement_5_source $path_file_supplement_5_product
cp $path_file_supplement_6_source $path_file_supplement_6_product
cp $path_file_supplement_7_source $path_file_supplement_7_product



#
