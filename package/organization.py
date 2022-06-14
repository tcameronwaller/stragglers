
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
    paths["bipolar_organization"] = os.path.join(
        path_dock, "bipolar_organization"
    )
    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["bipolar_organization"])
    # Initialize directories.
    utility.create_directories(
        path=paths["bipolar_organization"]
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
    path_table_phenotypes = os.path.join(
        path_dock, "bipolar_assembly",
        "table_phenotypes.pickle"
    )
    # Read information from file.
    table_phenotypes = pandas.read_pickle(
        path_table_phenotypes
    )
    # Compile and return information.
    return {
        "table_phenotypes": table_phenotypes,
    }


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





    # Collect information.
    pail_write = dict()
    pail_write["bipolar_organization"] = dict()
    pail_write["bipolar_organization"]["table_phenotypes"] = table
    # Write product information to file.
    write_product(
        pail_write=pail_write,
        paths=paths,
    )

    pass





#
