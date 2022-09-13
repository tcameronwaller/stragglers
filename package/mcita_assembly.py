
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
import uk_biobank.organization as ukb_org
import promiscuity.utility as utility
import promiscuity.polygenic_score as pgs
#import promiscuity.plot as plot

###############################################################################
# Functionality

# TODO: 1. I need an "access" script for the CITA data set components
# TODO: 1.1. collect and organize paths to A. phenotypes B. polygenic scores



# TODO: Read in the Polygenic Scores (PGS) like I did for the Bipolar Biobank.
# TODO:
#    # Read and organize tables of polygenic scores.
#    tables_polygenic_scores = pgs.drive_read_organize_tables_polygenic_scores(
#        table_scores_parameters=source["table_scores_parameters"],
#        filter_inclusion=True,
#        report=True,
#    )




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
    paths["mcita_assembly"] = os.path.join(path_dock, "mcita_assembly")
    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["mcita_assembly"])
    # Initialize directories.
    utility.create_directories(
        path=paths["mcita_assembly"]
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
        "polygenic_scores", "table_mayo_cita.tsv"
    )
    path_table_phenotypes = os.path.join(
        path_dock, "access", "mayo_cita",
        "table_phenotype.csv"
    )
    path_table_identifiers_case = os.path.join(
        path_dock, "access", "mayo_cita",
        "table_identifier_case.txt"
    )
    path_table_identifiers_control = os.path.join(
        path_dock, "access", "mayo_cita",
        "table_identifier_control.csv"
    )

    # Read information from file.
    table_parameter_scores = pgs.read_source_collection_polygenic_scores(
        path_table=path_table_parameter_scores,
        report=report,
    )

    table_phenotypes = pandas.read_csv(
        path_table_phenotypes,
        sep=",", # ","; "\t"; "\s+"; "\s+|\t+|\s+\t+|\t+\s+"
        header=0,
        dtype="string",
    )
    table_phenotypes.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )

    table_identifiers_case = pandas.read_csv(
        path_table_identifiers_case,
        sep="\s+", # ","; "\t"; "\s+"; "\s+|\t+|\s+\t+|\t+\s+"
        header=0,
        dtype="string",
    )
    table_identifiers_case.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )

    table_identifiers_control = pandas.read_csv(
        path_table_identifiers_control,
        sep=",", # ","; "\t"; "\s+"; "\s+|\t+|\s+\t+|\t+\s+"
        header=0,
        dtype="string",
    )
    table_identifiers_control.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )

    # Collect and return information.
    return {
        "table_parameter_scores": table_parameter_scores,
        "table_phenotypes": table_phenotypes,
        "table_identifiers_case": table_identifiers_case,
        "table_identifiers_control": table_identifiers_control,
    }


##########
# Organize separate tables before merge


def organize_table_column_identifier(
    column_source=None,
    column_product=None,
    table=None,
    report=None,
):
    """
    Organizes table of information about phenotypes.

    arguments:
        column_source (str): name of original column for identifier
        column_product (str): name of novel column to which to copy identifier
        table (object): Pandas data frame of information about phenotypes
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of information about phenotypes

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Convert all identifiers to type string.
    table[column_source] = table[column_source].astype("string")
    # Replace any empty identifier strings with missing values.
    table[column_source].replace(
        "",
        numpy.nan,
        inplace=True,
    )
    # Remove any records with missing identifiers.
    table.dropna(
        axis="index", # drop rows with missing values in columns
        how="any",
        subset=[column_source,],
        inplace=True,
    )
    # Convert identifiers to type string.
    # Copy identifiers.
    table[column_source] = table[column_source].astype("string")
    table[column_product] = table[column_source].astype("string").copy(
        deep=True,
    )
    # Remove columns.
    #table.drop(
    #    labels=["bib_id"],
    #    axis="columns",
    #    inplace=True
    #)
    # Return information.
    return table


def reduce_table_columns(
    columns_keep=None,
    table=None,
    report=None,
):
    """
    Simplifies or filters the columns in a table.

    arguments:
        columns_keep (list<str>): names of columns to keep in table
        table (object): Pandas data frame of information about phenotypes
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of information about phenotypes

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Filter columns in table.
    table = table.loc[
        :, table.columns.isin(columns_keep)
    ]
    # Return information.
    return table


def simplify_translate_table_columns_organize_identifier(
    columns_keep=None,
    columns_translations=None,
    columns_copy=None,
    identifier_source=None,
    identifier_product=None,
    table=None,
    report=None,
):
    """
    Organizes table of information about phenotypes.

    arguments:
        columns_keep (list<str>): names of columns to keep in table
        columns_translations (dict<str>): translation names for columns
        columns_copy (dict<str>): names of columns to copy
        identifier_source (str): name of original column for identifier
        identifier_product (str): name of novel column to which to copy identifier
        table (object): Pandas data frame of information about phenotypes
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of information about phenotypes

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Reduce, rename, and copy columns.
    table = reduce_table_columns(
        columns_keep=columns_keep,
        table=table,
        report=report,
    )
    table.rename(
        columns=columns_translations,
        inplace=True,
    )
    for column_new in columns_copy.keys():
        table[column_new] = table[columns_copy[column_new]].copy(deep=True)
    # Organize table index identifier.
    table = organize_table_column_identifier(
        column_source=identifier_source,
        column_product=identifier_product,
        table=table,
        report=report,
    )
    # Return information.
    return table



##########
# Genotype identifiers


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
        pail_write=pail_write["mcita_assembly"],
        path_directory=paths["mcita_assembly"],
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
    # Organize identifiers in tables before merges.
    table_phenotypes = organize_table_column_identifier(
        column_source="clinic_or_btogid",
        column_product="identifier_phenotype",
        table=source["table_phenotypes"],
        report=True,
    )

    # Convert identifiers to type string.
    # Copy identifiers.
    table_identifiers_case = (
        simplify_translate_table_columns_organize_identifier(
            columns_keep=["Sample.ID", "ClinicNum",],
            columns_translations={},
            columns_copy={
                "identifier_genotype_case": "Sample.ID"
            },
            identifier_source="ClinicNum",
            identifier_product="identifier_phenotype",
            table=source["table_identifiers_case"],
            report=True,
    ))
    # "RLIMS_Id", "External_Participant_Id", "genotype_dna_sampleid", "DNA_sampleid", "plasma_sampleid", "Sample_ID", "SubjectIden"
    table_identifiers_control = (
        simplify_translate_table_columns_organize_identifier(
            columns_keep=["control_btogid", "genotype_dna_sampleid",],
            columns_translations={},
            columns_copy={
                "identifier_genotype_control": "genotype_dna_sampleid",
            },
            identifier_source="control_btogid",
            identifier_product="identifier_phenotype",
            table=source["table_identifiers_control"],
            report=True,
    ))

    print("phenotype table before merges...")
    print(table_phenotypes)
    print("identifiers for cases...")
    print(table_identifiers_case)
    print("...")
    print("...")
    print("...")
    print("identifiers for controls...")
    print(table_identifiers_control)

    # Merge with phenotype variables the genotype identifiers for cases.
    table = utility.merge_columns_two_tables(
        identifier_first="identifier_phenotype",
        identifier_second="identifier_phenotype",
        table_first=table_phenotypes,
        table_second=table_identifiers_case,
        report=True,
    )
    # Merge with phenotype variables the genotype identifiers for controls.
    table = utility.merge_columns_two_tables(
        identifier_first="identifier_phenotype",
        identifier_second="identifier_phenotype",
        table_first=table,
        table_second=table_identifiers_control,
        report=True,
    )
    # Remove unnecessary columns from transformations on tables.
    table.drop(
        labels=["index_x", "index_y",],
        axis="columns",
        inplace=True
    )
    # Select relevant columns for phenotype variables.
    table = reduce_table_columns(
        columns_keep=[
            "clinic_or_btogid", "Sample.ID", "ClinicNum",
            "control_btogid", "genotype_dna_sampleid",
            "identifier_phenotype",
            "identifier_genotype_case", "identifier_genotype_control",
            "gender", "male", "age", "bmi", "race", "menopause", "bc", "case",
            "Albumin", "e1", "e2", "fsh", "lh", "Progesterone", "shbg",
            "Testosterone", "e2_", "e1_", "shbg_", "fsh_", "progest", "testost",
            "luteinizing", "albumin_",
            "log_e1", "log_e2", "log_fsh",
        ],
        table=table,
        report=True,
    )

    # Combine and organize genotype identifiers for cases and controls.
    table["identifier_genotype_case"] = (
        table["identifier_genotype_case"].astype("string")
    )
    table["identifier_genotype_control"] = (
        table["identifier_genotype_control"].astype("string")
    )
    table["identifier_genotype_case"].replace(
        numpy.nan,
        "",
        inplace=True,
    )
    table["identifier_genotype_control"].replace(
        numpy.nan,
        "",
        inplace=True,
    )
    # Determine consensus combination of identifiers for genotypes.
    table["identifier_genotype"] = table.apply(
        lambda row:
            utility.prioritize_combination_values_string(
                value_priority=row["identifier_genotype_control"],
                value_spare=row["identifier_genotype_case"],
            ),
        axis="columns", # apply function to each row
    )
    # Remove unnecessary columns from transformations on tables.
    table.drop(
        labels=["identifier_genotype_case", "identifier_genotype_control",],
        axis="columns",
        inplace=True
    )

    table = organize_table_column_identifier(
        column_source="identifier_genotype",
        column_product="identifier_genotype",
        table=table,
        report=True,
    )


    print("...")
    print("...")
    print("...")
    print("table after merge...")
    print(table)

    # Merge polygenic scores with information on phenotypes.
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
    pail_write["mcita_assembly"] = dict()
    pail_write["mcita_assembly"]["table_phenotypes"] = table
    # Write product information to file.
    write_product(
        pail_write=pail_write,
        paths=paths,
    )

    pass





#
