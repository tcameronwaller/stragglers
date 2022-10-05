
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
import promiscuity.plot as plot
import uk_biobank.stratification as ukb_strat

###############################################################################
# Functionality

# TODO: TCW; 5 October 2022
# TODO: preserve versatile functionality such as 'initialize_directories' and 'read_source' and 'write_product'


# TODO: TCW; 5 October 2022
# TODO: priorities
# 1. Test the new Transformations on Distribution Scale of variables AFTER stratification
# 1.1. stratify tables for cohorts
# 1.2. apply transformations
# 1.3. plot variables of interest within those cohorts



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
    paths["stragglers_scrap"] = os.path.join(path_dock, "stragglers_scrap")
    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["stragglers_scrap"])
    # Initialize directories.
    utility.create_directories(
        path=paths["stragglers_scrap"]
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
        path_dock, "organization_temp_copy",
        "table_phenotypes.pickle",
    )

    # Read information from file.
    table_phenotypes = pandas.read_pickle(
        path_table_phenotypes
    )
    # Compile and return information.
    return {
        "table_phenotypes": table_phenotypes,
        #"table_ukb_samples": table_ukb_samples,
    }



##########
# Write

# TODO: Follow pattern in 'ukbiobank.description'


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

    # Stratify records within separate tables for cohorts.
    records_cohorts = (
        ukb_strat.stratify_phenotype_cohorts_set_sex_age_menopause(
            table=source["table_phenotypes"],
    ))
    # Apply Distribution Scale Transformations to variables of interest in each
    # cohort.
    records_cohorts = (
        ukb_strat.stratify_phenotype_cohorts_set_sex_age_menopause(
            records_cohorts=records_cohorts,
    ))



    # Collect information.
    pail_write = dict()
    pail_write["stragglers_scrap"] = dict()
    pail_write["stragglers_scrap"]["tables_cohorts"] = tables_cohorts
    # Write product information to file.
    write_product(
        pail_write=pail_write,
        paths=paths,
    )

    pass





#
