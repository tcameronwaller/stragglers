
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
    paths["assembly_bipolar"] = os.path.join(path_dock, "assembly_bipolar")
    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["assembly_bipolar"])
    # Initialize directories.
    utility.create_directories(
        path=paths["assembly_bipolar"]
    )
    # Return information.
    return paths


##########
# Read


# TODO: TCW; 10 June 2022
# TODO: I want this to recognize ".LDpred2.auto" files and read them with a for loop
# TODO: also simplify the name of the column for the score...

def read_source(
    path_dock=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    Notice that Pandas does not accommodate missing values within series of
    integer variable types.

    The UK Biobank "eid" designates unique persons in the cohort.
    The UK Biobank "IID" matches persons to their genotype information.

    arguments:
        path_dock (str): path to dock directory for source and product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (object): source information

    """

    # Specify directories and files.
    path_table_phenotypes = os.path.join(
        path_dock, "access", "mayo_bipolar_phenotypes",
        "220513_BP_phenotypes.csv"
    )

    path_table_scores_steroid_globulin_female = os.path.join(
        path_dock, "access", "mayo_bipolar_polygenic_scores",
        "steroid_globulin_imputation_log_female_joint_1",
        "steroid_globulin_imputation_log_female_joint_1.LDpred2.auto"
    )
    path_table_scores_steroid_globulin_male = os.path.join(
        path_dock, "access", "mayo_bipolar_polygenic_scores",
        "steroid_globulin_imputation_log_male_joint_1",
        "steroid_globulin_imputation_log_male_joint_1.LDpred2.auto"
    )
    path_table_scores_testosterone_female = os.path.join(
        path_dock, "access", "mayo_bipolar_polygenic_scores",
        "testosterone_imputation_log_female_joint_1",
        "testosterone_imputation_log_female_joint_1.LDpred2.auto"
    )
    path_table_scores_testosterone_male = os.path.join(
        path_dock, "access", "mayo_bipolar_polygenic_scores",
        "testosterone_imputation_log_male_joint_1",
        "testosterone_imputation_log_male_joint_1.LDpred2.auto"
    )
    # Read information from file.
    table_phenotypes = pandas.read_csv(
        path_table_phenotypes,
        sep=",",
        header=0,
        dtype="string",
    )
    table_scores_steroid_globulin_female = pandas.read_csv(
        path_table_scores_steroid_globulin_female,
        sep="\t", # "\s+"; "\s+|\t+|\s+\t+|\t+\s+"
        header=0,
        dtype={
            "FID": "string",
            "IID": "string",
            "steroid_globulin_imputation_log_female_joint_1_auto": "float32",
        },
    )
    table_scores_steroid_globulin_male = pandas.read_csv(
        path_table_scores_steroid_globulin_male,
        sep="\t", # "\s+"; "\s+|\t+|\s+\t+|\t+\s+"
        header=0,
        dtype={
            "FID": "string",
            "IID": "string",
            "steroid_globulin_imputation_log_male_joint_1_auto": "float32",
        },
    )
    table_scores_testosterone_female = pandas.read_csv(
        path_table_scores_testosterone_female,
        sep="\t", # "\s+"
        header=0,
        dtype={
            "FID": "string",
            "IID": "string",
            "testosterone_imputation_log_female_joint_1_auto": "float32",
        },
    )
    table_scores_testosterone_male = pandas.read_csv(
        path_table_scores_testosterone_male,
        sep="\t", # "\s+"
        header=0,
        dtype={
            "FID": "string",
            "IID": "string",
            "testosterone_imputation_log_male_joint_1_auto": "float32",
        },
    )

    # Compile and return information.
    return {
        "table_phenotypes": table_phenotypes,
        "table_scores_steroid_globulin_female": table_scores_steroid_globulin_female,
        "table_scores_steroid_globulin_male": table_scores_steroid_globulin_male,
        "table_scores_testosterone_female": table_scores_testosterone_female,
        "table_scores_testosterone_male": table_scores_testosterone_male,
    }


##########
# Merge polygenic scores to phenotypes


def merge_polygenic_scores_to_phenotypes(
    table_phenotypes=None,
    table_scores_steroid_globulin_female=None,
    table_scores_steroid_globulin_male=None,
    table_scores_testosterone_female=None,
    table_scores_testosterone_male=None,
    report=None,
):
    """
    Removes irrelevant columns and rows from data.

    arguments:
        table_phenotypes (object): Pandas data frame of information about
            phenotype variables
        table_scores_... (object): Pandas data frame of polygenic scores
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of information about phenotype variables

    """

    # Copy data.
    table = table_phenotypes.copy(deep=True)

    # Translate names of columns.
    #table = table.add_prefix("import_")
    translations = dict()
    translations["steroid_globulin_imputation_log_female_joint_1_auto"] = (
        "steroid_globulin_female"
    )
    translations["steroid_globulin_imputation_log_male_joint_1_auto"] = (
        "steroid_globulin_male"
    )
    translations["testosterone_imputation_log_female_joint_1_auto"] = (
        "testosterone_female"
    )
    translations["testosterone_imputation_log_male_joint_1_auto"] = (
        "testosterone_male"
    )
    table_scores_steroid_globulin_female.rename(
        columns=translations,
        inplace=True,
    )
    table_scores_steroid_globulin_male.rename(
        columns=translations,
        inplace=True,
    )
    table_scores_testosterone_female.rename(
        columns=translations,
        inplace=True,
    )
    table_scores_testosterone_male.rename(
        columns=translations,
        inplace=True,
    )

    # Organize table indices.
    table.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    table.dropna(
        axis="index",
        how="any",
        subset=["bib_id"],
        inplace=True,
    )
    table["bib_id"].astype("string")
    table.set_index(
        "bib_id",
        append=False,
        drop=True,
        inplace=True
    )
    # Organize table indices.
    table_scores_steroid_globulin_female.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    table_scores_steroid_globulin_male.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    table_scores_testosterone_female.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    table_scores_testosterone_male.reset_index(
        level=None,
        inplace=True,
        drop=False,
    )
    table_scores_steroid_globulin_female.dropna(
        axis="index",
        how="any",
        subset=["IID"],
        inplace=True,
    )
    table_scores_steroid_globulin_male.dropna(
        axis="index",
        how="any",
        subset=["IID"],
        inplace=True,
    )
    table_scores_testosterone_female.dropna(
        axis="index",
        how="any",
        subset=["IID"],
        inplace=True,
    )
    table_scores_testosterone_male.dropna(
        axis="index",
        how="any",
        subset=["IID"],
        inplace=True,
    )
    table_scores_steroid_globulin_female["IID"].astype("string")
    table_scores_steroid_globulin_male["IID"].astype("string")
    table_scores_testosterone_female["IID"].astype("string")
    table_scores_testosterone_male["IID"].astype("string")
    table_scores_steroid_globulin_female.set_index(
        "IID",
        append=False,
        drop=True,
        inplace=True
    )
    table_scores_steroid_globulin_male.set_index(
        "IID",
        append=False,
        drop=True,
        inplace=True
    )
    table_scores_testosterone_female.set_index(
        "IID",
        append=False,
        drop=True,
        inplace=True
    )
    table_scores_testosterone_male.set_index(
        "IID",
        append=False,
        drop=True,
        inplace=True
    )
    table_scores_steroid_globulin_female.drop(
        labels=["FID"],
        axis="columns",
        inplace=True
    )
    table_scores_steroid_globulin_male.drop(
        labels=["FID"],
        axis="columns",
        inplace=True
    )
    table_scores_testosterone_female.drop(
        labels=["FID"],
        axis="columns",
        inplace=True
    )
    table_scores_testosterone_male.drop(
        labels=["FID"],
        axis="columns",
        inplace=True
    )
    # Merge data tables using database-style join.
    # Alternative is to use DataFrame.join().
    table = pandas.merge(
        table, # left table
        table_scores_steroid_globulin_female, # right table
        left_on="bib_id",
        right_on="IID",
        left_index=False,
        right_index=False,
        how="left", # keep only keys from left table
        suffixes=("_main", "_score"),
    )
    table = pandas.merge(
        table, # left table
        table_scores_steroid_globulin_male, # right table
        left_on="bib_id",
        right_on="IID",
        left_index=False,
        right_index=False,
        how="left", # keep only keys from left table
        suffixes=("_main", "_score"),
    )
    table = pandas.merge(
        table, # left table
        table_scores_testosterone_female, # right table
        left_on="bib_id",
        right_on="IID",
        left_index=False,
        right_index=False,
        how="left", # keep only keys from left table
        suffixes=("_main", "_score"),
    )
    table = pandas.merge(
        table, # left table
        table_scores_testosterone_male, # right table
        left_on="bib_id",
        right_on="IID",
        left_index=False,
        right_index=False,
        how="left", # keep only keys from left table
        suffixes=("_main", "_score"),
    )
    # Return information.
    return table



##########
# Remove irrelevant field columns


def parse_field_instance_columns_to_keep(
    instances_raw=None,
):
    """
    Parse the field's instances to keep.

    arguments:
        instances_raw (str): raw string of field's instances to keep

    raises:

    returns:
        (list<str>): field's instances to keep

    """

    instances_simple = instances_raw.replace("\'", "")
    instances_simple = instances_simple.replace("\"", "")
    instances = list()
    for instance_raw in instances_simple.split(";"):
        instance = str(instance_raw).strip()
        if len(instance) > 0:
            instances.append(str(instance))
            pass
        pass
    return instances


def parse_ukbiobank_raw_column_field_instance_array(
    name_column=None,
):
    """
    Parse the information in a raw column name.

    Format: [field]-[instance].[array index]
    Examples: "20002-3.33", "20003-3.47"

    arguments:
        name_column (str): raw name of column in data table from UK Biobank

    raises:

    returns:
        (dict<str>): information from column name

    """

    record = dict()
    if (
        ("-" in name_column) and
        ("." in name_column)
    ):
        record["field"] = str(name_column.split("-")[0].strip())
        record["instance"] = str(
            name_column.split("-")[1].split(".")[0].strip()
        )
        record["array_index"] = str(
            name_column.split("-")[1].split(".")[1].strip()
        )
    else:
        record["field"] = ""
        record["instance"] = ""
        record["array_index"] = ""
    return record


def determine_keep_column_field_instance(
    column=None,
    table_ukbiobank_variables=None,
    report=None,
):
    """
    Determines whether to keep a column for UK Biobank field instance.

    arguments:
        column (str): column name from UK Biobank accessions
        table_ukbiobank_variables (object): Pandas data frame of information
            about UK Biobank phenotype variables
        report (bool): whether to print reports

    raises:

    returns:
        (bool): whether to keep the column

    """

    # Parse information in the column name.
    record = parse_ukbiobank_raw_column_field_instance_array(
        name_column=column,
    )
    # Determine whether column matches format of an original accession
    # data-field.
    if (len(str(record["field"])) > 0):
        # Initialize flag.
        keep = False
        # Organize information.
        #array_collection = table_ukbiobank_variables.at[
        #    int(record["field"]), "array_collection",
        #]
        instances_keep_raw = table_ukbiobank_variables.at[
            int(record["field"]), "instances_keep",
        ]
        # Determine which instances to keep.
        if not pandas.isna(instances_keep_raw):
            # Organize field instances to keep.
            instances_keep = parse_field_instance_columns_to_keep(
                instances_raw=instances_keep_raw,
            )
            # Determine whether to keep column on basis of instance.
            if str(record["instance"]) in instances_keep:
                keep = True
                # Report.
                if (
                    report and
                    (str(record["field"]) in ["31", "50", "22009"])
                ):
                    utility.print_terminal_partition(level=4)
                    print("report from: determine_keep_column_field_instance()")
                    print("keep column: " + str(column))
                pass
            pass
        pass
    else:
        # Column is not an original accession data-field.
        # Keep the column.
        keep = True
    # Return information.
    return keep


def determine_ukbiobank_field_instance_columns_keep(
    columns_accession=None,
    table_ukbiobank_variables=None,
    extra_names=None,
    report=None,
):
    """
    Organizes column names for variable fields and instances.

    arguments:
        columns_accession (list<str>): column names from UK Biobank accessions
        table_ukbiobank_variables (object): Pandas data frame of information
            about UK Biobank phenotype variables
        extra_names (list<str>): extra names to include
        report (bool): whether to print reports

    raises:

    returns:
        (list<str>): column names for variable fields and instances

    """

    # Copy data.
    table_ukbiobank_variables = table_ukbiobank_variables.copy(deep=True)
    # Organize information.
    table_ukbiobank_variables = table_ukbiobank_variables.loc[
        :, table_ukbiobank_variables.columns.isin([
            "field", "array_collection", "instances_keep"
        ])
    ]
    table_ukbiobank_variables["field"].astype("string")
    table_ukbiobank_variables.set_index(
        "field",
        drop=True,
        inplace=True,
    )
    # Iterate on actual column names from the accession table.
    # Determine whether to keep column.
    columns_keep = list()
    columns_keep.extend(extra_names)
    for column in columns_accession:
        column = str(column)
        # Determine whether to keep column.
        keep = determine_keep_column_field_instance(
            column=column,
            table_ukbiobank_variables=table_ukbiobank_variables,
            report=False,
        )
        if keep:
            columns_keep.append(column)
            pass
        pass
    # Return information.
    return columns_keep


def remove_table_irrelevant_field_instance_columns(
    table_ukbiobank_variables=None,
    columns_accession=None,
    table_ukb_41826=None,
    table_ukb_43878=None,
    table_ukb_47488=None,
    report=None,
):
    """
    Removes irrelevant columns and rows from data.

    arguments:
        table_ukbiobank_variables (object): Pandas data frame of information
            about UK Biobank phenotype variables
        columns_accession (list<str>): column names from UK Biobank accessions
        table_ukb_41826 (object): Pandas data frame of variables from UK
            Biobank phenotype data accession 41826
        table_ukb_43878 (object): Pandas data frame of variables from UK
            Biobank phenotype data accession 43878
        table_ukb_47488 (object): Pandas data frame of variables from UK
            Biobank phenotype data accession 47488
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of Pandas data frames after removal of
            irrelevant columns and rows

    """

    # Copy data.
    table_variables = table_ukbiobank_variables.copy(deep=True)
    table_ukb_41826 = table_ukb_41826.copy(deep=True)
    table_ukb_43878 = table_ukb_43878.copy(deep=True)
    table_ukb_47488 = table_ukb_47488.copy(deep=True)
    # Determine which columns to keep for UK Biobank fields and instances.
    columns_keep = determine_ukbiobank_field_instance_columns_keep(
        columns_accession=columns_accession,
        table_ukbiobank_variables=table_variables,
        extra_names=["IID", "eid"],
        report=report,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("columns to keep: " + str(len(columns_keep)))
        #print(columns_keep[0:50])
        print(columns_keep)
        utility.print_terminal_partition(level=3)
        print("...before pruning...")
        print("table_ukb_41826 shape: " + str(table_ukb_41826.shape))
        print("table_ukb_43878 shape: " + str(table_ukb_43878.shape))
        print("table_ukb_47488 shape: " + str(table_ukb_47488.shape))
        utility.print_terminal_partition(level=4)
    # Remove all irrelevant columns.
    table_ukb_41826 = table_ukb_41826.loc[
        :, table_ukb_41826.columns.isin(columns_keep)
    ]
    table_ukb_43878 = table_ukb_43878.loc[
        :, table_ukb_43878.columns.isin(columns_keep)
    ]
    table_ukb_47488 = table_ukb_47488.loc[
        :, table_ukb_47488.columns.isin(columns_keep)
    ]
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("...after pruning...")
        print("table_ukb_41826 shape: " + str(table_ukb_41826.shape))
        print("table_ukb_43878 shape: " + str(table_ukb_43878.shape))
        print("table_ukb_47488 shape: " + str(table_ukb_47488.shape))
        utility.print_terminal_partition(level=4)
    # Collect information.
    pail = dict()
    pail["table_ukb_41826"] = table_ukb_41826
    pail["table_ukb_43878"] = table_ukb_43878
    pail["table_ukb_47488"] = table_ukb_47488
    # Return information.
    return pail


##########
# Organization


def simplify_field_values_array_row(
    row=None,
    field=None,
    columns_field=None,
    delimiter=None,
    report=None,
):
    """
    Simplify field instances for array values.

    arguments:
        row (object): Pandas data frame row
        field (str): identifier of UK Biobank field
        columns_field (list<str>): names of columns that match the data-field
        delimiter (str): delimiter for string representation of array values
        report (bool): whether to print reports

    raises:

    returns:
        (str): text array with semicolon delimiters

    """

    # Copy Pandas series for row.
    row = row.copy(deep=True)
    # Convert Pandas series row to dictionary.
    record_row = row.to_dict()
    # Collect all non-missing values from the field's columns for instances.
    values = list()
    for column in columns_field:
        # Parse information in the column name.
        record_column = parse_ukbiobank_raw_column_field_instance_array(
            name_column=column,
        )
        # Determine whether column matches format of an original accession
        # data-field.
        if (
            (column in record_row.keys()) and
            (len(str(record_column["field"])) > 0) and
            (str(field) == str(record_column["field"]))
        ):
            value = record_row[column]
            if (
                (not pandas.isna(value)) and
                (str(value) != "<NA>")
            ):
                # Collect all values regardless of whether they are unique.
                values.append(str(value))
                pass
            pass
        pass
    # Select unique values.
    values_unique = list(set(values)) # unique
    # Combine values with text delimiter.
    if len(values_unique) > 0:
        text_array = delimiter.join(values_unique)
    else:
        text_array = ""
    # Report.
    if report:
        utility.print_terminal_partition(level=4)
        print("matching fields to " + str(field) + " :")
        print(columns_field)
        print("unique values: " + str(len(values_unique)))
    # Return information.
    return text_array


def simplify_field_values_array_columns(
    table_ukbiobank_variables=None,
    table_ukb_raw=None,
    delimiter=None,
    report=None,
):
    """
    Simplify field instances for array values.

    arguments:
        table_ukbiobank_variables (object): Pandas data frame of information
            about UK Biobank phenotype variables
        table_ukb_raw (object): Pandas data frame of information from UK
            Biobank
        delimiter (str): delimiter for string representation of array values
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of information from UK Biobank

    """

    # Copy and de-fragment information in table.
    table_ukb = table_ukb_raw.copy(deep=True)
    table_variables = table_ukbiobank_variables.copy(deep=True)
    # Organize information.
    table_variables["field"].astype("string")
    table_variables = table_variables.loc[
        :, table_variables.columns.isin(["field", "type", "array_collection"])
    ]
    table_variables["array_collection"] = table_variables.apply(
        lambda row:
            str(row["array_collection"]).strip().lower(),
        axis="columns", # apply across rows
    )
    table_variables_array = table_variables.loc[
        (
            ~pandas.isna(table_variables["array_collection"]) &
            (table_variables["array_collection"] == "yes")
        ), :
    ]
    fields_array = table_variables_array["field"].to_list()
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("fields for arrays to simplify: " + str(len(fields_array)))
        print(fields_array)
    # Iterate on UK Biobank fields with array instances.
    for field in fields_array:
        # Find columns that match the data-field.
        columns_field = list()
        for column in table_ukb.columns.to_list():
            # Parse information in the column name.
            record = parse_ukbiobank_raw_column_field_instance_array(
                name_column=column,
            )
            # Determine whether column matches format of an original accession
            # data-field.
            # Determine whether the column matches the data-field.
            if (
                (len(str(record["field"])) > 0) and
                (str(field) == str(record["field"]))
            ):
                columns_field.append(column)
                pass
            pass
        # Create new column with text array of all non-missing values from the
        # field's original columns.
        column_new = str(str(field) + "_array")
        table_ukb[column_new] = table_ukb.apply(
            lambda row:
                simplify_field_values_array_row(
                    row=row,
                    field=field,
                    columns_field=columns_field,
                    delimiter=delimiter,
                    report=False,
                ),
            axis="columns", # apply across rows
        )
        # Report.
        if report:
            utility.print_terminal_partition(level=3)
            print("dropping columns for field: " + str(field))
            print(columns_field)
        # Drop columns.
        table_ukb.drop(
            labels=columns_field,
            axis="columns",
            inplace=True
        )
        pass
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("...after simplification of array fields...")
        print("table_ukb shape: " + str(table_ukb.shape))
        utility.print_terminal_partition(level=4)
    # Return information.
    return table_ukb


def merge_table_variables_identifiers(
    table_identifier_pairs=None,
    table_ukb_41826=None,
    table_ukb_43878=None,
    table_ukb_47488=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    arguments:
        table_identifier_pairs (object): Pandas data frame of associations
            between "IID" and "eid"
        table_ukb_41826 (object): Pandas data frame of variables from UK Biobank
            phenotype accession 41826
        table_ukb_43878 (object): Pandas data frame of variables from UK Biobank
            phenotype accession 43878
        table_ukb_47488 (object): Pandas data frame of variables from UK Biobank
            phenotype accession 47488
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(table_identifier_pairs)
        utility.print_terminal_partition(level=2)
        print(table_ukb_41826)
        utility.print_terminal_partition(level=2)
        print(table_ukb_43878)
        utility.print_terminal_partition(level=2)
        print(table_ukb_47488)
    # Remove rows with null values of merge identifier.
    table_identifier_pairs.dropna(
        axis="index",
        how="any",
        subset=["eid"],
        inplace=True,
    )
    table_ukb_41826.dropna(
        axis="index",
        how="any",
        subset=["eid"],
        inplace=True,
    )
    table_ukb_43878.dropna(
        axis="index",
        how="any",
        subset=["eid"],
        inplace=True,
    )
    table_ukb_47488.dropna(
        axis="index",
        how="any",
        subset=["eid"],
        inplace=True,
    )
    # Organize data.
    table_identifier_pairs.astype("string")
    table_identifier_pairs.set_index(
        "eid",
        drop=True,
        inplace=True,
    )
    table_ukb_41826["eid"].astype("string")
    table_ukb_41826.set_index(
        "eid",
        drop=True,
        inplace=True,
    )
    table_ukb_43878["eid"].astype("string")
    table_ukb_43878.set_index(
        "eid",
        drop=True,
        inplace=True,
    )
    table_ukb_47488["eid"].astype("string")
    table_ukb_47488.set_index(
        "eid",
        drop=True,
        inplace=True,
    )

    # Merge data tables using database-style join.
    # Alternative is to use DataFrame.join().
    table_merge = table_identifier_pairs.merge(
        table_ukb_41826,
        how="outer",
        left_on="eid",
        right_on="eid",
        suffixes=("_merge", "_41826"),
    )
    table_merge = table_merge.merge(
        table_ukb_43878,
        how="outer",
        left_on="eid",
        right_on="eid",
        suffixes=("_merge", "_43878"),
    )
    table_merge = table_merge.merge(
        table_ukb_47488,
        how="outer",
        left_on="eid",
        right_on="eid",
        suffixes=("_merge", "_47488"),
    )
    # Remove excess columns.

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print(table_merge)
    # Return information.
    return table_merge


def exclude_persons_ukbiobank_consent(
    exclusion_identifiers=None,
    table=None,
    report=None,
):
    """
    Removes entries with specific identifiers from data.

    arguments:
        exclusion_identifiers (list<str>): identifiers of persons who withdrew
            consent from UK Biobank
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of phenotype variables across UK Biobank
            cohort

    """

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("Initial table dimensions: " + str(table.shape))
        utility.print_terminal_partition(level=4)
        print("Exclusion of persons: " + str(len(exclusion_identifiers)))
    # Copy data.
    table = table.copy(deep=True)
    # Filter data entries.
    table_exclusion = table.loc[
        ~table.index.isin(exclusion_identifiers), :
    ]
    # Report.
    if report:
        utility.print_terminal_partition(level=4)
        print("Final data dimensions: " + str(table_exclusion.shape))
    # Return information.
    return table_exclusion


def drop_null_records_all_variables(
    table=None,
    columns_any=None,
    report=None,
):
    """
    Drops rows with null or missing values across all columns.

    arguments:
        table (object): Pandas data frame of information about UK Biobank
            phenotype variables
        columns_any (str): column names for which to drop rows with any missing
            values
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of information about UK Biobank phenotype
            variables

    """

    # Copy data.
    table = table.copy(deep=True)
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("...before dropping null records...")
        print("table shape: " + str(table.shape))
        print(table)
        utility.print_terminal_partition(level=4)
    # Remove rows with null values across all columns.
    table.dropna(
        axis="index",
        how="all",
        inplace=True,
    )
    table.dropna(
        axis="index",
        how="any",
        subset=columns_any,
        inplace=True,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("...after dropping null records...")
        print("table shape: " + str(table.shape))
        print(table)
        utility.print_terminal_partition(level=4)
    # Return information.
    return table


##########
# Write


def write_product_raw(
    information=None,
    path_parent=None,
):
    """
    Writes product information to file.

    arguments:
        information (object): information to write to file
        path_parent (str): path to parent directory
            raises:

    returns:

    """

    # Specify directories and files.
    path_table_ukb_41826_text = os.path.join(
        path_parent, "table_ukb_41826.tsv"
    )
    # Write information to file.
    information["table_ukb_41826"].to_csv(
        path_or_buf=path_table_ukb_41826_text,
        sep="\t",
        header=True,
        index=True,
    )
    pass


def write_product_inspection(
    information=None,
    path_parent=None,
):
    """
    Writes product information to file.

    arguments:
        information (object): information to write to file
        path_parent (str): path to parent directory
            raises:

    returns:

    """

    # Specify directories and files.
    path_table_text = os.path.join(
        path_parent, "table_phenotypes.tsv"
    )
    # Write information to file.
    information["table_phenotypes"].to_csv(
        path_or_buf=path_table_text,
        sep="\t",
        header=True,
        index=True,
    )
    pass


def write_product_assembly(
    information=None,
    path_parent=None,
):
    """
    Writes product information to file.

    arguments:
        information (object): information to write to file
        path_parent (str): path to parent directory
            raises:

    returns:

    """

    # Specify directories and files.
    path_table_phenotypes = os.path.join(
        path_parent, "table_phenotypes.pickle"
    )
    path_table_phenotypes_text = os.path.join(
        path_parent, "table_phenotypes.tsv"
    )
    path_table_kinship_pairs = os.path.join(
        path_parent, "table_kinship_pairs.pickle"
    )
    path_table_kinship_pairs_text = os.path.join(
        path_parent, "table_kinship_pairs.tsv"
    )
    # Write information to file.
    information["table_phenotypes"].to_pickle(
        path_table_phenotypes
    )
    information["table_phenotypes"].to_csv(
        path_or_buf=path_table_phenotypes_text,
        sep="\t",
        header=True,
        index=True,
    )
    information["table_kinship_pairs"].to_pickle(
        path_table_kinship_pairs
    )
    information["table_kinship_pairs"].to_csv(
        path_or_buf=path_table_kinship_pairs_text,
        sep="\t",
        header=True,
        index=False,
    )
    pass


def write_product(
    information=None,
    paths=None,
):
    """
    Writes product information to file.

    arguments:
        information (object): information to write to file
        paths (dict<str>): collection of paths to directories for procedure's
            files

    raises:

    returns:

    """

    write_product_raw(
        information=information["raw"],
        path_parent=paths["raw"],
    )
    write_product_inspection(
        information=information["inspection"],
        path_parent=paths["inspection"],
    )
    write_product_assembly(
        information=information["assembly"],
        path_parent=paths["assembly"],
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
    # Exclusion identifiers are "eid".
    source = read_source(
        path_dock=path_dock,
        report=True,
    )

    # Mege polygenic scores with information on phenotypes.
    table = merge_polygenic_scores_to_phenotypes(
        table_phenotypes=source["table_phenotypes"],
        table_scores_steroid_globulin_female=source["table_scores_steroid_globulin_female"],
        table_scores_steroid_globulin_male=source["table_scores_steroid_globulin_male"],
        table_scores_testosterone_female=source["table_scores_testosterone_female"],
        table_scores_testosterone_male=source["table_scores_testosterone_male"],
        report=True,
    )

    if False:
        # Remove data columns for irrelevant variable instances.
        prune = remove_table_irrelevant_field_instance_columns(
            table_ukbiobank_variables=source["table_ukbiobank_variables"],
            columns_accession=source["columns_accession"],
            table_ukb_41826=source["table_ukb_41826"],
            table_ukb_43878=source["table_ukb_43878"],
            table_ukb_47488=source["table_ukb_47488"],
            report=True,
        )

        # Simplify UK Biobank fields with multiple instances.
        # Reduce these multiple field instance columns to arrays.
        table_ukb_41826_simple = simplify_field_values_array_columns(
            table_ukbiobank_variables=source["table_ukbiobank_variables"],
            table_ukb_raw=prune["table_ukb_41826"],
            delimiter=";",
            report=True,
        )

        # Merge tables.
        table_merge = merge_table_variables_identifiers(
            table_identifier_pairs=source["table_identifier_pairs"],
            table_ukb_41826=table_ukb_41826_simple,
            table_ukb_43878=prune["table_ukb_43878"],
            table_ukb_47488=prune["table_ukb_47488"],
            report=True,
        )

        # Exclude persons who withdrew consent from the UK Biobank.
        table_exclusion = exclude_persons_ukbiobank_consent(
            exclusion_identifiers=source["exclusion_identifiers"],
            table=table_merge,
            report=True,
        )
        # Drop any records (persons) with null values across all variables.
        table_valid = drop_null_records_all_variables(
            table=table_exclusion,
            columns_any=["31-0.0"], # "31-0.0": sex
            report=True,
        )

        utility.print_terminal_partition(level=2)
        print("table_kinship_pairs")
        print(source["table_kinship_pairs"])

        # Write out raw tables for inspection.
        # Collect information.
        information = dict()
        information["raw"] = dict()
        information["inspection"] = dict()
        information["assembly"] = dict()
        information["raw"]["table_ukb_41826"] = (
            source["table_ukb_41826"].iloc[0:10000, :]
        )
        information["inspection"]["table_phenotypes"] = (
            table_valid.iloc[0:10000, :]
        )
        information["assembly"]["table_phenotypes"] = table_valid
        information["assembly"]["table_kinship_pairs"] = (
            source["table_kinship_pairs"]
        )
        # Write product information to file.
        write_product(
            paths=paths,
            information=information
        )
    pass





#
