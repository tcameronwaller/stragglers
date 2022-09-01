#!/bin/bash

#chmod u+x script.sh

################################################################################
# Notes:
# first draft of script: TCW; 27 July 2022
# last revision of script: TCW; 27 July 2022

################################################################################
# Organize paths to files and directories.

# Read private, local file path for main project directory.
cd ~/paths
path_directory_process=$(<"./process_sexy_alcohol.txt")
path_directory_storage=$(<"./storage_waller_metabolism.txt")
path_directory_case_control_data=$(<"./mayo_cita_case_control_data.txt")
path_directory_case_data=$(<"./mayo_cita_case_data.txt")

# Define paths.
path_directory_dock="${path_directory_process}/dock"
path_directory_persistence="${path_directory_storage}/phenotypes/mayo_clinic_cita_access_2022-07-27"
path_directory_access="${path_directory_dock}/access/mayo_cita"
path_file_phenotype_source="${path_directory_storage}/phenotypes/mayo_clinic_cita_access_2022-07-27/AUDhormone_analysis.csv"
path_file_identifier_control_source="${path_directory_case_control_data}/Clin_geneIIDs_04152022.csv" # file date: 15 April 2022
#path_file_identifier_control_source="${path_directory_case_control_data}/winham_Karpyak_case_control_match_QN_data_12022021.csv"
path_file_identifier_case_source="${path_directory_case_control_data}/CITAcaseIID_to_SubjIden.txt" # file date: 28 January 2022
#path_file_identifier_case_source="${path_directory_case_data}/rlim_id_map_2021-11-29_1548.csv"
path_file_identifier_control_persistence="${path_directory_persistence}/Clin_geneIIDs_04152022.csv"
path_file_identifier_case_persistence="${path_directory_persistence}/CITAcaseIID_to_SubjIden.txt"
path_file_phenotype_product="${path_directory_access}/table_phenotype.csv"
path_file_identifier_control_product="${path_directory_access}/table_identifier_control.csv"
path_file_identifier_case_product="${path_directory_access}/table_identifier_case.txt"

# Echo each command to console.
set -x

# Remove previous versions.
rm -r $path_directory_access

# Create directory.
mkdir -p $path_directory_access

################################################################################
# Copy files to persistence directory.

if false; then
  cp $path_file_identifier_control_source $path_file_identifier_control_persistence
  cp $path_file_identifier_case_source $path_file_identifier_case_persistence
fi

################################################################################
# Copy files to processing directory.

cp $path_file_phenotype_source $path_file_phenotype_product
cp $path_file_identifier_control_source $path_file_identifier_control_product
cp $path_file_identifier_case_source $path_file_identifier_case_product
