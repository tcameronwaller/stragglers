
"""

This module contains functions for general organization of the raw extraction
data from the UK Biobank accessions.

Combine information from multiple UK Biobank accessions.
Integrate identifier matchings.
Exclude individual cases (persons) who subsequently withdrew consent.
It is also necessary to handle the somewhat awkward format of array fields.

"""

###############################################################################
# Notes

###############################################################################
# Installation and importation

# Import modules from specific path without having to install a general package
# I would have to figure out how to pass a path variable...
# https://stackoverflow.com/questions/67631/how-to-import-a-module-given-the-full-path


# Standard

import sys
#print(sys.path)
import os
import math
import statistics
import pickle
import copy
import random
import itertools
import time

# Relevant

import numpy
import pandas
pandas.options.mode.chained_assignment = None # default = "warn"

import scipy.stats
import scipy.linalg
import statsmodels.multivariate.pca

# Custom
import promiscuity.utility as utility
import promiscuity.polygenic_score as pgs
import stragglers.mcita_assembly

###############################################################################
# Functionality


##########
# Initialization


def initialize_directories(
    restore=None,
    path_dock=None,
):
    """
    Initialize directories for procedure's product files.

    arguments:
        restore (bool): whether to remove previous versions of data
        path_dock (str): path to dock directory for source and product
            directories and files

    raises:

    returns:
        (dict<str>): collection of paths to directories for procedure's files

    """

    # Collect paths.
    paths = dict()
    # Define paths to directories.
    paths["dock"] = path_dock
    paths["bipolar_assembly"] = os.path.join(path_dock, "bipolar_assembly")
    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["bipolar_assembly"])
    # Initialize directories.
    utility.create_directories(
        path=paths["bipolar_assembly"]
    )
    # Return information.
    return paths


##########
# Read


def read_source(
    path_dock=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    arguments:
        path_dock (str): path to dock directory for source and product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_table_identifiers = os.path.join(
        path_dock, "access", "mayo_bipolar_phenotypes",
        "210421_id_matching_gwas.csv"
    )
    path_table_phenotypes = os.path.join(
        path_dock, "access", "mayo_bipolar_phenotypes",
        "220513_BP_phenotypes.csv"
    )
    path_table_genetic_sex_case = os.path.join(
        path_dock, "access", "mayo_bipolar_phenotypes",
        "MERGED.maf0.01.dosR20.8.noDups.fam"
    )
    path_table_parameter_scores = os.path.join(
        path_dock, "parameters", "bipolar_biobank",
        "polygenic_scores", "table_bipolar_biobank.tsv"
    )
    # Read information from file.
    table_parameter_scores = pgs.read_source_collection_polygenic_scores(
        path_table=path_table_parameter_scores,
        report=report,
    )
    table_identifiers = pandas.read_csv(
        path_table_identifiers,
        sep=",",
        header=0,
        #dtype="string",
        dtype={
            "bib_id": "string",
            "gwas1_sampleid": "string", # identifier of individual's genotype
            "gwas2_sampleid": "string", # identifier of individual's genotype
        },
    )
    table_identifiers.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_phenotypes = pandas.read_csv(
        path_table_phenotypes,
        sep=",",
        header=0,
        dtype="string",
    )
    table_phenotypes.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    # https://www.cog-genomics.org/plink/2.0/formats#fam
    table_genetic_sex_case = pandas.read_csv(
        path_table_genetic_sex_case,
        sep="\s+", # "\t"; "\s+"; "\s+|\t+|\s+\t+|\t+\s+"
        header=None,
        names=[
            "FID", "IID", "father", "mother",
            "sex_genotype_raw", "bipolar_disorder_genotype_raw"
        ],
        dtype={
            "FID": "string",
            "IID": "string", # identifier of individual's genotype
            "father": "string",
            "mother": "string",
            "sex_genotype_raw": "string", # 1: male; 2: female; 0: unknown
            "bipolar_disorder_genotype_raw": "string", # 1: control; 2: case;
        },
    )
    table_genetic_sex_case.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    # Compile and return information.
    return {
        "table_identifiers": table_identifiers,
        "table_phenotypes": table_phenotypes,
        "table_genetic_sex_case": table_genetic_sex_case,
        "table_parameter_scores": table_parameter_scores,
    }


##########
# Organize separate tables before merge


def prioritize_genotype_identifiers(
    genotype_identifier_priority=None,
    genotype_identifier_spare=None,
):
    """
    Determines the identifier for a priority genotype record.

    arguments:
        genotype_identifier_priority (str): identifier for a priority genotype
        genotype_identifier_spare (str): identifier for a spare genotype

    raises:

    returns:
        (str): identifier for priority genotype

    """

    # Determine identifier of priority genotype.
    if (
        (not pandas.isna(genotype_identifier_priority)) and
        (len(str(genotype_identifier_priority)) > 0)
    ):
        # Identifier for priority genotype record is not missing.
        # Priority genotype record comes from a primary, priority set or
        # batch of genotypes.
        genotype_priority = str(
            copy.deepcopy(genotype_identifier_priority)
        )
    elif (
        (not pandas.isna(genotype_identifier_spare)) and
        (len(str(genotype_identifier_spare)) > 0)
    ):
        # Identifier for spare genotype record is not missing.
        # Spare genotype record comes from a secondary, not priority set or
        # batch of genotypes.
        genotype_priority = str(
            copy.deepcopy(genotype_identifier_spare)
        )
        pass
    else:
        # There is not a genotype record available to match the phenotype
        # record.
        genotype_priority = ""
    # Return information.
    return genotype_priority


def organize_table_phenotype_genotype_identifiers(
    table=None,
    report=None,
):
    """
    Organizes table of identifiers for phenotype and genotype records.

    arguments:
        table (object): Pandas data frame of identifiers for matching phenotype
            and genotype records
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of identifiers for matching phenotype and
            genotype records

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Convert all identifiers to type string.
    table["bib_id"] = table["bib_id"].astype("string")
    table["gwas1_sampleid"] = table["gwas1_sampleid"].astype("string")
    table["gwas2_sampleid"] = table["gwas2_sampleid"].astype("string")
    # Prioritize identifiers from "GWAS1" set of genotypes.
    table["identifier_genotype"] = table.apply(
        lambda row:
            prioritize_genotype_identifiers(
                genotype_identifier_priority=row["gwas1_sampleid"],
                genotype_identifier_spare=row["gwas2_sampleid"],
            ),
        axis="columns", # apply function to each row
    )
    # Replace any empty identifier strings with missing values.
    table["bib_id"].replace(
        "",
        numpy.nan,
        inplace=True,
    )
    table["identifier_genotype"].replace(
        "",
        numpy.nan,
        inplace=True,
    )
    # Remove any records with missing identifiers.
    table.dropna(
        axis="index", # drop rows with missing values in columns
        how="all",
        subset=["bib_id", "identifier_genotype"],
        inplace=True,
    )
    # Convert identifiers to type string.
    table["bib_id"] = table["bib_id"].astype("string")
    table["identifier_phenotype"] = table["bib_id"].astype("string")
    table["identifier_genotype"] = table["identifier_genotype"].astype("string")
    # Return information.
    return table


def organize_table_genetic_sex_case(
    table=None,
    report=None,
):
    """
    Organizes table of information about phenotypes.

    arguments:
        table (object): Pandas data frame of information about phenotypes
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of information about phenotypes

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Convert all identifiers to type string.
    table["IID"] = table["IID"].astype("string")
    # Replace any empty identifier strings with missing values.
    table["IID"].replace(
        "",
        numpy.nan,
        inplace=True,
    )
    # Remove any records with missing identifiers.
    table.dropna(
        axis="index", # drop rows with missing values in columns
        how="any",
        subset=["IID",],
        inplace=True,
    )
    # Convert identifiers to type string.
    table["IID"] = table["IID"].astype("string")
    table["identifier_genotype"] = table["IID"].astype("string")
    # Remove the column for the genotype family identifier ("FID").
    # Remove the columns for identifiers of parents.
    # These identifiers are redundant and unnecessary.
    table.drop(
        labels=["FID", "IID", "father", "mother"],
        axis="columns",
        inplace=True
    )
    # Return information.
    return table


##########
# Merge information on phenotypes


# TODO: split up this function... some of this function is specific to the Bipolar Biobank... joining phenotypes, identifiers, etc...


def merge_polygenic_scores_to_phenotypes(
    table_identifiers=None,
    table_phenotypes=None,
    table_genetic_sex_case=None,
    tables_scores=None,
    report=None,
):
    """
    Merge information from multiple source tables.

    Most controls that accompany the Bipolar Biobank do not have phenotype
    records and do not have a "bib_id".
    Only a few controls from the Mexico and Chile cohorts do have "bib_id" and
    phenotype records.

    Many cases in the Bipolar Biobank do not have genotypes.

    Sample identifiers from genotype files ("IID") are the most inclusive
    identifiers when analyses prioritize genotypes.

    arguments:
        table_identifiers (object): Pandas data frame of identifiers for
            matching phenotype and genotype records
        table_phenotypes (object): Pandas data frame of information about
            phenotype variables
        table_genetic_sex_case (object): Pandas data frame of information about
            genetic sex and case-control status from file in PLINK2 ".fam"
            format
        tables_scores (list<object>): collection of Pandas data frames of
            polygenic scores
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of information about phenotype variables

    """

    # 1. Introduce phenotype records for Bipolar Disorder cases (mostly) to the
    # table of genetic sex and case status.
    # The table of genetic sex and case status is the most inclusive of cases
    # and controls with genotypes.

    # 1.1. Introduce identifiers for phenotype records ("bib_id") to the table
    # of genetic sex and case status.
    # Copy information in table.
    table_genetic_sex_case = table_genetic_sex_case.copy(deep=True)
    table_identifiers = table_identifiers.copy(deep=True)
    # Organize tables' indices.
    table_genetic_sex_case.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_genetic_sex_case["identifier_genotype"] = (
        table_genetic_sex_case["identifier_genotype"].astype("string")
    )
    table_genetic_sex_case.set_index(
        "identifier_genotype",
        append=False,
        drop=True, # move regular column to index; remove original column
        inplace=True
    )
    table_identifiers.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_identifiers["identifier_genotype"] = (
        table_identifiers["identifier_genotype"].astype("string")
    )
    table_identifiers.set_index(
        "identifier_genotype",
        append=False,
        drop=True, # move regular column to index; remove original column
        inplace=True
    )
    # Merge data tables using database-style join.
    # Alternative is to use DataFrame.join().
    table = pandas.merge(
        table_genetic_sex_case, # left table
        table_identifiers, # right table
        left_on=None, # "identifier_genotype",
        right_on=None, # "identifier_genotype",
        left_index=True,
        right_index=True,
        how="outer", # keep union of keys from both tables
        #suffixes=("_main", "_identifiers"), # deprecated?
    )

    # 1.2. Introduce phenotype records for Bipolar Disorder cases (mostly) to
    # the main table of genetic sex and case status.
    # Copy information in table.
    table_phenotypes = table_phenotypes.copy(deep=True)
    # Organize tables' indices.
    table.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    table["identifier_phenotype"] = (
        table["identifier_phenotype"].astype("string")
    )
    table.set_index(
        "identifier_phenotype",
        append=False,
        drop=True, # move regular column to index; remove original column
        inplace=True
    )
    table_phenotypes.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_phenotypes["identifier_phenotype"] = (
        table_phenotypes["identifier_phenotype"].astype("string")
    )
    table_phenotypes.set_index(
        "identifier_phenotype",
        append=False,
        drop=True, # move regular column to index; remove original column
        inplace=True
    )
    # Merge data tables using database-style join.
    # Alternative is to use DataFrame.join().
    table = pandas.merge(
        table, # left table
        table_phenotypes, # right table
        left_on=None, # "identifier_phenotype",
        right_on=None, # "identifier_phenotype",
        left_index=True,
        right_index=True,
        how="outer", # keep union of keys from both tables
        #suffixes=("_main", "_identifiers"), # deprecated?
    )

    # 3. Introduce polygenic scores.
    # Organize main table's index.
    table.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    table["identifier_genotype"] = (
        table["identifier_genotype"].astype("string")
    )
    table.set_index(
        "identifier_genotype",
        append=False,
        drop=True, # move regular column to index; remove original column
        inplace=True
    )
    # Iterate on tables for polygenic scores.
    for table_score in tables_scores:
        # Organize score table's index.
        table_score.reset_index(
            level=None,
            inplace=True,
            drop=True, # remove index; do not move to regular columns
        )
        table_score["identifier_genotype"] = (
            table_score["identifier_genotype"].astype("string")
        )
        table_score.set_index(
            "identifier_genotype",
            append=False,
            drop=True, # move regular column to index; remove original column
            inplace=True
        )
        # Merge data tables using database-style join.
        # Alternative is to use DataFrame.join().
        table = pandas.merge(
            table, # left table
            table_score, # right table
            left_on=None, # "identifier_genotype",
            right_on=None, # "identifier_genotype",
            left_index=True,
            right_index=True,
            how="outer", # keep union of keys from both tables
            #suffixes=("_main", "_identifiers"), # deprecated?
        )
        pass
    # Organize table's index.
    table.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    #table.set_index(
    #    "identifier_genotype",
    #    append=False,
    #    drop=True, # move regular column to index; remove original column
    #    inplace=True
    #)
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        print("merge_polygenic_scores_to_phenotypes()")
        utility.print_terminal_partition(level=3)
        print(table)
        print("columns")
        print(table.columns.to_list())
        pass
    # Return information.
    return table


##########
# Organize phenotype variables

#    # Convert column types to float.
#    columns_type = [
#        "2784-0.0", "2794-0.0", "2804-0.0",
#        "2814-0.0", "3536-0.0", "3546-0.0",
#    ]
#    table = utility.convert_table_columns_variables_types_float(
#        columns=columns_type,
#        table=table,
#    )


##########
# Write


def write_product_bipolar_assembly(
    pail_write=None,
    path_directory=None,
):
    """
    Writes product information to file.

    arguments:
        pail_write (dict): collection of information to write to file
        path_directory (str): path to parent directory

    raises:

    returns:

    """

    # Specify directories and files.
    path_table_phenotypes = os.path.join(
        path_directory, "table_phenotypes.pickle"
    )
    path_table_phenotypes_text = os.path.join(
        path_directory, "table_phenotypes.tsv"
    )
    # Write information to file.
    pail_write["table_phenotypes"].to_pickle(
        path_table_phenotypes
    )
    pail_write["table_phenotypes"].to_csv(
        path_or_buf=path_table_phenotypes_text,
        sep="\t",
        header=True,
        index=False, # include index in table
    )
    pass


def write_product(
    pail_write=None,
    paths=None,
):
    """
    Writes product information to file.

    arguments:
        pail_write (dict): collection of information to write to file
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    # Organization procedure main information.
    write_product_bipolar_assembly(
        pail_write=pail_write["bipolar_assembly"],
        path_directory=paths["bipolar_assembly"],
    )
    pass


###############################################################################
# Procedure

# TODO: TCW; 13 June 2022
# TODO: 1. I need to read in the ancestry and/or ethnicity ("European" or "Hispanic")


def execute_procedure(
    path_dock=None,
):
    """
    Function to execute module's main behavior.

    arguments:
        path_dock (str): path to dock directory for source and product
            directories and files

    raises:

    returns:

    """

    # Report version.
    utility.print_terminal_partition(level=1)
    print(path_dock)
    print("version check: 1")
    # Pause procedure.
    time.sleep(5.0)

    # Initialize directories.
    paths = initialize_directories(
        restore=True,
        path_dock=path_dock,
    )
    # Read source information from file.
    source = read_source(
        path_dock=path_dock,
        report=True,
    )
    # Read and organize tables of polygenic scores.
    tables_polygenic_scores = pgs.drive_read_organize_tables_polygenic_scores(
        table_parameter_scores=source["table_parameter_scores"],
        filter_inclusion=True,
        report=True,
    )

    # Organize table of identifiers for phenotypes and genotypes.
    table_identifiers = organize_table_phenotype_genotype_identifiers(
        table=source["table_identifiers"],
        report=True,
    )
    # Organize table of phenotypes.
    table_phenotypes = (
        stragglers.mcita_assembly.organize_table_column_identifier(
            column_source="bib_id",
            column_target="identifier_phenotype",
            table=source["table_phenotypes"],
            report=True,
    ))
    # Organize table of genetic sex (from PLINK2 file in ".fam" format).
    # https://www.cog-genomics.org/plink/2.0/formats#fam
    table_genetic_sex_case = organize_table_genetic_sex_case(
        table=source["table_genetic_sex_case"],
        report=True,
    )
    # Merge together information about phenotypes.
    table_phenotypes_merge = merge_phenotype_variables_identifiers(
        table_identifiers=table_identifiers,
        table_phenotypes=table_phenotypes,
        table_genetic_sex_case=table_genetic_sex_case,
        report=True,
    )

    # Merge polygenic scores with information on phenotypes.
    table = utility.merge_tables_supplements_to_main(
        identifier_main="identifier_genotype",
        identifier_supplement="identifier_genotype",
        table_main=table_phenotypes_merge,
        tables_supplements=tables_polygenic_scores,
        report=True,
    )
    # TODO: drop redundant columns

    # Collect information.
    pail_write = dict()
    pail_write["bipolar_assembly"] = dict()
    pail_write["bipolar_assembly"]["table_phenotypes"] = table
    # Write product information to file.
    write_product(
        pail_write=pail_write,
        paths=paths,
    )

    pass





#
