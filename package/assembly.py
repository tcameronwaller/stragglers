
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
#import promiscuity.plot as plot

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
    path_table_scores_parameters = os.path.join(
        path_dock, "parameters", "bipolar_biobank",
        "polygenic_scores", "table_polygenic_scores.tsv"
    )
    # Read information from file.
    table_identifiers = pandas.read_csv(
        path_table_identifiers,
        sep=",",
        header=0,
        dtype="string",
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
    table_scores_parameters = pandas.read_csv(
        path_table_scores_parameters,
        sep="\t", # "\t"; "\s+"; "\s+|\t+|\s+\t+|\t+\s+"
        header=0,
        dtype={
            "inclusion": "int32",
            "name_file_path_directory_parent": "string",
            "path_directory": "string",
            "name_file": "string",
            "method": "string",
            "name_column_score_old": "string",
            "name_column_score_new": "string",
        },
    )
    table_scores_parameters.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    # Compile and return information.
    return {
        "table_identifiers": table_identifiers,
        "table_phenotypes": table_phenotypes,
        "table_scores_parameters": table_scores_parameters,
    }


##########
# Read and organize tables of polygenic scores


def read_source_table_polygenic_score_ldpred2(
    name_file_directory_parent=None,
    path_directory=None,
    name_file=None,
    name_column_score=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    arguments:
        name_file_directory_parent (str): name of file with private path to
            parent directory
        path_directory (str): path beyond parent directory to the file
        name_file (str): name of file
        name_column_score (str): name of column in table for polygenic scores
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of polygenic scores

    """

    # Read private path to parent directory.
    path_home_paths = os.path.expanduser("~/paths")
    path_file_directory_parent = os.path.join(
        path_home_paths, name_file_directory_parent
    )
    path_directory_parent = utility.read_file_text(
        path_file=path_file_directory_parent
    ).rstrip("\n") # remove new line character from string
    # Specify directories and files.
    path_table = os.path.join(
        path_directory_parent, path_directory, name_file
    )
    # Read information from file.
    table = pandas.read_csv(
        path_table,
        sep="\s+", # "\t"; "\s+"; "\s+|\t+|\s+\t+|\t+\s+"
        header=0,
        dtype={
            "FID": "string",
            "IID": "string",
            name_column_score: "float32",
        },
    )
    table.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    # Compile and return information.
    return table


def organize_table_polygenic_score_ldpred2(
    table=None,
    name_score_old=None,
    name_score_new=None,
    report=None,
):
    """
    Organizes table of polygenic scores.

    arguments:
        table (object): Pandas data frame of information about polygenic scores
        name_score_old (str): old, original name for polygenic score variable
        name_score_new (str): new, novel name for polygenic score variable
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
    # Translate names of columns.
    #table = table.add_prefix("import_")
    translations = dict()
    translations[name_score_old] = name_score_new
    table.rename(
        columns=translations,
        inplace=True,
    )
    # Convert scores to type float.
    #table[name_score_new] = table[name_score_new].astype("float32")
    table[name_score_new] = pandas.to_numeric(
        table[name_score_new],
        errors="coerce", # force any invalid values to missing or null
        downcast="float",
    )
    # Remove the column for the genotype family identifier ("FID").
    # This identifier is redundant and unnecessary.
    table.drop(
        labels=["FID"],
        axis="columns",
        inplace=True
    )
    # Return information.
    return table


def drive_read_organize_tables_polygenic_scores(
    table_scores_parameters=None,
    filter_inclusion=None,
    report=None,
):
    """
    Drives functions to read and organize source information from file.

    arguments:
        table_scores_parameters (object): Pandas data frame of parameters for
           reading and organizing separate tables for polygenic scores
        filter_inclusion (bool): whether to filter records in polygenic scores
            parameter table by logical binary "inclusion" variable
        report (bool): whether to print reports

    raises:

    returns:
        (list<object>): collection of Pandas data frames of polygenic scores

    """

    # Filter collection of polygenic scores by "inclusion" flag.
    if filter_inclusion:
        table_scores_parameters = table_scores_parameters.loc[
            (
                (table_scores_parameters["inclusion"] == 1)
            ), :
        ]
        pass
    # Extract records from table.
    records = table_scores_parameters.to_dict(
        orient="records",
    )
    # Collect tables of polygenic scores.
    tables_polygenic_scores = list()
    # Iterate on records.
    for record in records:
        # Extract information from record.
        name_file_directory_parent = record["name_file_path_directory_parent"]
        path_directory = record["path_directory"]
        name_file = record["name_file"]
        method = record["method"]
        name_column_score_old = record["name_column_score_old"]
        name_column_score_new = record["name_column_score_new"]
        # Determine appropriate functions for format of polygenic scores.
        if (method == "LDPred2"):
            # Read file.
            table_raw = read_source_table_polygenic_score_ldpred2(
                name_file_directory_parent=name_file_directory_parent,
                path_directory=path_directory,
                name_file=name_file,
                name_column_score=name_column_score_old,
                report=report,
            )
            # Organize table.
            table = organize_table_polygenic_score_ldpred2(
                table=table_raw,
                name_score_old=name_column_score_old,
                name_score_new=name_column_score_new,
                report=report,
            )
            # Collect table.
            tables_polygenic_scores.append(table)
        elif (method == "PRS-CS"):
            print("still need to implement read and organize for PRS-CS...")
            pass
        # Report.
        if report:
            utility.print_terminal_partition(level=2)
            print("report: ")
            print("drive_read_organize_tables_polygenic_scores()")
            utility.print_terminal_partition(level=3)
            print(table)
            print("columns")
            print(table.columns.to_list())
            pass
    # Return information.
    return tables_polygenic_scores



##########
# Organize separate tables before merge


def prioritize_genotype_identifiers(
    phenotype_identifier=None,
    genotype_identifier_priority=None,
    genotype_identifier_spare=None,
):
    """
    Determines the identifier for a priority genotype record that matches the
    identifier for a phenotype record.

    arguments:
        phenotype_identifer (str): identifier for a phenotype record
        genotype_identifier_priority (str): identifier for a priority genotype
        genotype_identifier_spare (str): identifier for a spare genotype

    raises:

    returns:
        (str): identifier for priority genotype

    """

    # Determine the priority genotype record that matches the phenotype record.
    if (
        (not pandas.isna(phenotype_identifier)) and
        (len(str(phenotype_identifier)) > 0)
    ):
        # Identifier for phenotype record is not missing.
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
        pass
    else:
        # The is not a phenotype record available.
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
                phenotype_identifier=row["bib_id"],
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
        how="any",
        subset=["bib_id", "identifier_genotype"],
        inplace=True,
    )
    # Convert identifiers to type string.
    table["bib_id"] = table["bib_id"].astype("string")
    table["identifier_genotype"] = table["identifier_genotype"].astype("string")
    # Return information.
    return table


def organize_table_phenotypes(
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
    table["bib_id"] = table["bib_id"].astype("string")
    # Replace any empty identifier strings with missing values.
    table["bib_id"].replace(
        "",
        numpy.nan,
        inplace=True,
    )
    # Remove any records with missing identifiers.
    table.dropna(
        axis="index", # drop rows with missing values in columns
        how="any",
        subset=["bib_id",],
        inplace=True,
    )
    # Convert identifiers to type string.
    table["bib_id"] = table["bib_id"].astype("string")
    # Return information.
    return table


##########
# Merge polygenic scores to phenotypes


def merge_polygenic_scores_to_phenotypes(
    table_identifiers=None,
    table_phenotypes=None,
    tables_scores=None,
    report=None,
):
    """
    Removes irrelevant columns and rows from data.

    arguments:
        table_identifiers (object): Pandas data frame of identifiers for
            matching phenotype and genotype records
        table_phenotypes (object): Pandas data frame of information about
            phenotype variables
        tables_scores (list<object>): collection of Pandas data frames of
            polygenic scores
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of information about phenotype variables

    """

    # 1. Introduce identifiers of genotype records to table of phenotype
    # records.
    # Copy information in table.
    table_identifiers = table_identifiers.copy(deep=True)
    table_phenotypes = table_phenotypes.copy(deep=True)
    # Organize tables' indices.
    table_identifiers.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_identifiers["bib_id"] = table_identifiers["bib_id"].astype("string")
    table_identifiers.set_index(
        "bib_id",
        append=False,
        drop=True, # move regular column to index; remove original column
        inplace=True
    )
    table_phenotypes.reset_index(
        level=None,
        inplace=True,
        drop=True, # remove index; do not move to regular columns
    )
    table_phenotypes["bib_id"] = table_phenotypes["bib_id"].astype("string")
    table_phenotypes.set_index(
        "bib_id",
        append=False,
        drop=True, # move regular column to index; remove original column
        inplace=True
    )
    # Merge data tables using database-style join.
    # Alternative is to use DataFrame.join().
    table = pandas.merge(
        table_phenotypes, # left table
        table_identifiers, # right table
        left_on=None, # "bib_id",
        right_on=None, # "IID",
        left_index=True,
        right_index=True,
        how="left", # keep only keys from left table
        #suffixes=("_main", "_identifiers"), # deprecated?
    )

    # 2. Introduce polygenic scores.
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
        drop=False, # move regular column to index; remove original column
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
        table_score["IID"] = table_score["IID"].astype("string")
        table_score.set_index(
            "IID",
            append=False,
            drop=True, # move regular column to index; remove original column
            inplace=True
        )
        # Merge data tables using database-style join.
        # Alternative is to use DataFrame.join().
        table = pandas.merge(
            table, # left table
            table_score, # right table
            left_on=None, # "bib_id",
            right_on=None, # "IID",
            left_index=True,
            right_index=True,
            how="left", # keep only keys from left table
            #suffixes=("_main", "_score"), # deprecated?
        )
        pass
    # Organize table's index.
    table.reset_index(
        level=None,
        inplace=True,
        drop=False, # remove index; do not move to regular columns
    )
    table.set_index(
        "bib_id",
        append=False,
        drop=True, # move regular column to index; remove original column
        inplace=True
    )
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
        index=True,
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
# TODO: 1. I need to read in the genetic sex and ancestry and or ethnicity ("European" or "Hispanic")
# TODO: 2. translate Bipolar Disorder diagnosis type
# TODO: 3. clean up rapid cycling variable
# TODO: 4. organize chronotype variables?



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
    # Exclusion identifiers are "eid".
    source = read_source(
        path_dock=path_dock,
        report=True,
    )
    # Read and organize tables of polygenic scores.
    tables_polygenic_scores = drive_read_organize_tables_polygenic_scores(
        table_scores_parameters=source["table_scores_parameters"],
        filter_inclusion=True,
        report=True,
    )

    # Organize table of identifiers for phenotypes and genotypes.
    table_identifiers = organize_table_phenotype_genotype_identifiers(
        table=source["table_identifiers"],
        report=True,
    )
    # Organize table of phenotypes.
    table_phenotypes = organize_table_phenotypes(
        table=source["table_phenotypes"],
        report=True,
    )
    # Mege polygenic scores with information on phenotypes.
    table = merge_polygenic_scores_to_phenotypes(
        table_identifiers=table_identifiers,
        table_phenotypes=table_phenotypes,
        tables_scores=tables_polygenic_scores,
        report=True,
    )

    # Collect information.
    pail_write = dict()
    pail_write["bipolar_assembly"] = dict()
    pail_write["bipolar_assembly"]["table_phenotypes"] = table
    # Write product information to file.
    write_product(
        pail_write=pail_write,
        paths=paths,
    )

    if False:

        # Organize phenotype variables.

        # "bib_id": phenotype identifier
        # "gender": gender
        # "pt_age": age
        # "BMI": body mass index
        # "rc" rapid cycling encoded as a binary variable (derived from multiple categories)
        # "scid_dx": Bipolar Disorder type I or II
        # "database": name of source database for phenotype (clinical) records
        # "SITE": assessment center?

        # Organize table.
        # Select relevant columns from table.
        columns_selection = [
            "bib_id",
            "identifier_genotype",
            "gender",
            "pt_age",
            "BMI",
            "rc",
            "scid_dx",
            "database",
            "steroid_globulin_female",
            "steroid_globulin_male",
            "testosterone_female",
            "testosterone_male",
        ]
        table = table.loc[
            :, table.columns.isin(columns_selection)
        ]
        utility.print_terminal_partition(level=2)
        print("table after selection of columns")
        print(table)
        print("columns")
        print(table.columns.to_list())

        table = table[[*columns_selection]]
        # Define logical binary indicator variables for type of Bipolar Disorder
        # diagnosis.
        table = define_logical_binary_indicator_variables_bipolar_disorder_type(
            table=table,
            report=True,
        )

    pass





#
