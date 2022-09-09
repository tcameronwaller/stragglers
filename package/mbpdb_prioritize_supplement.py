
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
import stragglers.mcita_assembly as s_mcita_ass

###############################################################################
# Functionality

# TODO: TCW; 9 September 2022
# TODO: 1. read in the phenotype table and all of the 6 candidate supplement tables
# TODO: 2. determine whether "pt_age" and "BMI" are same variable names in all tables
# TODO: - - This might require opening each .tsv file in Excel
# TODO: 3. stratify tables by Bipolar Disorder cases and controls
# TODO: 4. summarize missingness of Age and BMI variables in cases and controls in each candidate supplement table



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
    paths["mbpdb_assembly"] = os.path.join(path_dock, "mbpdb_assembly")
    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["mbpdb_assembly"])
    # Initialize directories.
    utility.create_directories(
        path=paths["mbpdb_assembly"]
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
    path_table_parameter_scores = os.path.join(
        path_dock, "parameters", "stragglers",
        "polygenic_scores", "table_mayo_bpdb.tsv"
    )

    path_table_identifiers = os.path.join(
        path_dock, "access", "mayo_bpdb",
        "210421_id_matching_gwas.csv"
    )
    path_table_phenotypes = os.path.join(
        path_dock, "access", "mayo_bpdb",
        "220513_BP_phenotypes.csv"
    )
    path_table_genetic_sex_case = os.path.join(
        path_dock, "access", "mayo_bpdb",
        "MERGED.maf0.01.dosR20.8.noDups.fam"
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



##########
# Merge information on phenotypes


# TODO: split up this function... some of this function is specific to the Bipolar Biobank... joining phenotypes, identifiers, etc...

# TODO: TCW; 01 September 2022
# TODO: obsolete
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


def write_product_assembly(
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
    write_product_assembly(
        pail_write=pail_write["mbpdb_assembly"],
        path_directory=paths["mbpdb_assembly"],
    )
    pass


###############################################################################
# Procedure


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

    # Organize table of phenotypes.
    table_phenotypes = (
        s_mcita_ass.organize_table_column_identifier(
            column_source="bib_id",
            column_product="identifier_phenotype",
            table=source["table_phenotypes"],
            report=True,
    ))

    # Organize table of identifiers.
    # Determine consensus combination of identifiers for genotypes.
    # Prioritize identifiers from "GWAS1" set of genotypes.
    table_identifiers = (
        s_mcita_ass.simplify_translate_table_columns_organize_identifier(
            columns_keep=["bib_id", "gwas1_sampleid", "gwas2_sampleid",],
            columns_translations={},
            columns_copy={},
            identifier_source="bib_id",
            identifier_product="identifier_phenotype",
            table=source["table_identifiers"],
            report=True,
    ))
    table_identifiers["gwas_sampleid_consensus"] = table_identifiers.apply(
        lambda row:
            s_mcita_ass.prioritize_combination_values(
                value_priority=row["gwas1_sampleid"],
                value_spare=row["gwas2_sampleid"],
            ),
        axis="columns", # apply function to each row
    )
    table_identifiers = s_mcita_ass.organize_table_column_identifier(
            column_source="gwas_sampleid_consensus",
            column_product="identifier_genotype",
            table=table_identifiers,
            report=True,
    )

    # Organize table of genetic sex (from PLINK2 file in ".fam" format).
    # https://www.cog-genomics.org/plink/2.0/formats#fam
    table_genetic_sex_case = (
        s_mcita_ass.simplify_translate_table_columns_organize_identifier(
            columns_keep=[
                "IID", "sex_genotype_raw", "bipolar_disorder_genotype_raw"
            ],
            columns_translations={},
            columns_copy={},
            identifier_source="IID",
            identifier_product="identifier_genotype",
            table=source["table_genetic_sex_case"],
            report=True,
    ))

    # Merge table of phenotype variables with table of phenotype and genotype
    # identifiers.
    table = utility.merge_columns_two_tables(
        identifier_first="identifier_phenotype",
        identifier_second="identifier_phenotype",
        table_first=table_phenotypes,
        table_second=table_identifiers,
        report=True,
    )
    # Merge table of phenotype variables with table of genetic sex and case
    # status.
    table = s_mcita_ass.organize_table_column_identifier(
            column_source="identifier_genotype",
            column_product="identifier_genotype",
            table=table,
            report=True,
    )
    table = utility.merge_columns_two_tables(
        identifier_first="identifier_genotype",
        identifier_second="identifier_genotype",
        table_first=table,
        table_second=table_genetic_sex_case,
        report=True,
    )
    # Remove unnecessary columns from transformations on tables.
    table.drop(
        labels=["index_x", "index_y",],
        axis="columns",
        inplace=True
    )

    # Merge table of phenotype variables with tables of polygenic scores.
    table = utility.merge_columns_tables_supplements_to_main(
        identifier_main="identifier_genotype",
        identifier_supplement="identifier_genotype",
        table_main=table,
        tables_supplements=tables_polygenic_scores,
        report=True,
    )
    # Remove unnecessary columns from transformations on tables.
    table.drop(
        labels=["index_x", "index_y", "index",],
        axis="columns",
        inplace=True
    )

    # Report.
    print("...")
    print("...")
    print("...")
    print("table after merges with PGS...")
    print(table)
    utility.print_terminal_partition(level=3)
    print("table columns: " + str(int(table.shape[1])))
    print("table rows: " + str(int(table.shape[0])))
    print("columns")
    print(table.columns.to_list())


    # Collect information.
    pail_write = dict()
    pail_write["mbpdb_assembly"] = dict()
    pail_write["mbpdb_assembly"]["table_phenotypes"] = table
    # Write product information to file.
    write_product(
        pail_write=pail_write,
        paths=paths,
    )

    pass





#
