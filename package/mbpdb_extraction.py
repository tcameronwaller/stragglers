
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
import partner.utility as utility
import partner.extraction as pextr
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
    if True:
        paths["correlation"] = os.path.join(
            path_dock, "ldsc_gwas_disorders_tcw_2023-08-02", "6_gwas_correlation_ldsc",
        )
        paths["extraction_ldsc"] = os.path.join(
            path_dock, "ldsc_gwas_disorders_tcw_2023-08-02", "6_gwas_correlation_ldsc_extraction",
        )
    if False:
        paths["heritability"] = os.path.join(
            path_dock, "alcohol_genetics_tcw_2023-03-02", "gwas_heritability_observed_liability_scale_ldsc",
        )
        paths["correlation"] = os.path.join(
            path_dock, "alcohol_genetics_tcw_2023-03-02", "gwas_correlation_ldsc",
        )
        paths["extraction_ldsc"] = os.path.join(
            path_dock, "alcohol_genetics_tcw_2023-03-02", "extraction_ldsc",
        )
    if False:
        paths["heritability"] = os.path.join(
            path_dock, "hormone_genetics_tcw_2023-02-24", "gwas_heritability_ldsc",
        )
        paths["correlation"] = os.path.join(
            path_dock, "hormone_genetics_tcw_2023-02-24", "gwas_correlation_ldsc",
        )
        paths["extraction_ldsc"] = os.path.join(
            path_dock, "hormone_genetics_tcw_2023-02-24", "extraction_ldsc",
        )
    if False:
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
        utility.remove_directory(path=paths["extraction_ldsc"])
    # Initialize directories.
    utility.create_directories(
        path=paths["extraction_ldsc"]
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

    # Collection information.
    pail_write = dict()

    # Extract names of child directories within parent directory.
    names_directories = utility.extract_subdirectory_names(
        path=paths["correlation"]
    )
    names_directories_ldsc = list(filter(
        lambda name: (name != "logs"),
        names_directories
    ))
    print("--------------------")
    print(names_directories_ldsc)
    print("--------------------")

    # Write each table to file.
    for name_directory in names_directories_ldsc:

        path_directory = os.path.join(
            paths["correlation"], name_directory,
        )

        table_correlation = pextr.read_extract_from_all_ldsc_files_in_directory(
            path_directory=path_directory,
            file_name_pattern=".log",
            file_name_pattern_not=".....",
            analysis="correlation",
            report=True,
        )
        pail_write[str("table_" + name_directory)] = table_correlation



    # Extract information from reports of analyses in LDSC.
    if False:
        table_heritability = pextr.read_extract_from_all_ldsc_files_in_directory(
            path_directory=paths["heritability"],
            file_name_pattern=".log",
            file_name_pattern_not=".....",
            analysis="heritability",
            report=True,
        )
        pail_write["table_heritability"] = table_heritability
    if False:
        table_correlation = pextr.read_extract_from_all_ldsc_files_in_directory(
            path_directory=paths["correlation"],
            file_name_pattern=".log",
            file_name_pattern_not=".....",
            analysis="correlation",
            report=True,
        )
        pail_write["table_correlation"] = table_correlation

    # Write information to file.
    control_write_product(
        pail_write=pail_write,
        path_directory=paths["extraction_ldsc"],
    )



    pass





#
