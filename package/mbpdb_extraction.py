
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
import promiscuity.extraction as pextr
import stragglers.mcita_assembly as s_mcita_ass

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
    paths["heritability"] = os.path.join(
        path_dock, "bipolar_body", "gwas_heritability_ldsc",
    )
    paths["correlation"] = os.path.join(
        path_dock, "bipolar_body", "gwas_genetic_correlation_ldsc",
    )
    paths["extraction"] = os.path.join(
        path_dock, "bipolar_body", "extraction_ldsc",
    )
    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["extraction"])
    # Initialize directories.
    utility.create_directories(
        path=paths["extraction"]
    )
    # Return information.
    return paths



##########
# Write


def write_product_table(
    name=None,
    table=None,
    path_directory=None,
):
    """
    Writes product information in a Pandas data frame table to file.

    arguments:
        name (str): base name for file
        table (object): Pandas data frame table
        path_directory (str): path to directory in which to write file

    raises:

    returns:

    """

    # Specify directories and files.
    path_table_text = os.path.join(
        path_directory, str(name + ".tsv")
    )
    # Write information to file.
    table.to_csv(
        path_or_buf=path_table_text,
        sep="\t",
        header=True,
        index=False,
        na_rep="NA",
    )
    pass


def control_write_product(
    pail_write=None,
    path_directory=None,
):
    """
    Writes product information to file.

    arguments:
        pail_write (dict): collection of information to write to file
        path_parent (str): full path to directory in which to write files

    raises:

    returns:

    """

    # Initialize directories.
    utility.create_directories(
        path=path_directory,
    )
    # Write each table to file.
    for name in pail_write.keys():
        write_product_table(
            name=name,
            table=pail_write[name],
            path_directory=path_directory,
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

    # Extract information from reports of analyses in LDSC.
    table_heritability = pextr.read_extract_from_all_ldsc_files_in_directory(
        path_directory=paths["heritability"],
        file_name_pattern=".log",
        file_name_pattern_not=".....",
        analysis="heritability",
        report=True,
    )
    table_correlation = pextr.read_extract_from_all_ldsc_files_in_directory(
        path_directory=paths["correlation"],
        file_name_pattern=".log",
        file_name_pattern_not=".....",
        analysis="correlation",
        report=True,
    )

    # Write product information to file.
    pail_write = dict()
    pail_write["table_heritability"] = table_heritability
    pail_write["table_correlation"] = table_correlation
    control_write_product(
        pail_write=pail_write,
        path_directory=paths["extraction"],
    )



    pass





#
