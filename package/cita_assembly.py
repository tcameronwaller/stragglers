
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
    paths["cita_assembly"] = os.path.join(path_dock, "cita_assembly")
    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["cita_assembly"])
    # Initialize directories.
    utility.create_directories(
        path=paths["cita_assembly"]
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
        path_dock, "access", "mayo_cita_phenotypes",
        "AUDhormone_analysis.csv"
    )
    # Read information from file.
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
    # Return information.
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



# Main driver


def drive_estimate_bioavailable_free_estradiol_testosterone(
    table=None,
    report=None,
):
    """
    Drives the estimation of bioavailable and free estradiol and testosterone.

    arguments:
        table (object): Pandas data frame table of phenotype variables
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data-frame table of phenotype variables

    """

    # Copy information in table.
    table = table.copy(deep=True)

    ##########
    # Organize raw hormone variables.

    # Convert variable types to float.
    columns_float = list()
    columns_float.append("e1_")
    columns_float.append("e2_")
    columns_float.append("testost")
    columns_float.append("shbg_")
    columns_float.append("albumin_")
    table = utility.convert_table_columns_variables_types_float(
        columns=columns_float,
        table=table,
    )
    # Convert concentrations to units of moles per liter (mol/L) with adjustment
    # by specific factors for appropriate scale in analyses (and floats).
    factors_concentration = dict()
    factors_concentration["oestradiol"] = 1E12 # 1 pmol / L
    factors_concentration["oestradiol_free"] = 1E12 # 1 pmol / L
    factors_concentration["oestradiol_bioavailable"] = 1E12 # 1 pmol / L
    factors_concentration["testosterone"] = 1E12 # 1 pmol / L
    factors_concentration["testosterone_free"] = 1E12 # 1 pmol / L
    factors_concentration["testosterone_bioavailable"] = 1E12 # 1 pmol / L
    factors_concentration["steroid_globulin"] = 1E9 # 1 nmol / L
    factors_concentration["albumin"] = 1E6 # 1 umol / L
    table = convert_hormone_concentration_units_moles_per_liter(
        table=table,
        factors_concentration=factors_concentration,
    )

    ##########
    # Calculate estimates of bioavailable and free hormones.
    table = organize_calculation_estimate_bioavailable_free_hormones(
        factors_concentration=factors_concentration,
        table=table,
        report=report,
    )


    # Return information.
    return table


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

# TODO: TCW; 06 July 2022
# TODO: 1. convert measurement variables to float32
# TODO: 2. convert measurement variables to units of molarity (moles / liter)
# TODO: 3. call the functions from "uk_biobank.organization" to estimate bioavailable and free testosterone and estradiol

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

    print(source["table_phenotypes"])

    # Calculate estimates of bioavailable and free estradiol and testosterone.
    table = drive_estimate_bioavailable_free_estradiol_testosterone(
        table=source["table_phenotypes"],
        report=True,
    )



    pass





#
