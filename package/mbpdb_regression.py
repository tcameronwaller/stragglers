"""
Organize regression analyses on data from the Bipolar Biobank.

Author:

    T. Cameron Waller
    tcameronwaller@gmail.com
    Rochester, Minnesota 55904
    United States of America

License:

    This file is part of bipolar_biobank
    (https://github.com/tcameronwaller/bipolar_biobank/).

    bipolar_biobank supports analyses on data from the Bipolar Biobank.
    Copyright (C) 2021 Thomas Cameron Waller

    bipolar_biobank is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option)
    any later version.

    UK_Biobank is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
    more details.

    You should have received a copy of the GNU General Public License along
    with UK_Biobank. If not, see <http://www.gnu.org/licenses/>.
"""

###############################################################################
# Notes

###############################################################################
# Installation and importation

# Standard

import sys
#print(sys.path)
import os
#import shutil
#import csv
import math
import statistics
import pickle
import copy
import random
import itertools
import textwrap
import time

# Relevant

import numpy
import scipy.stats
import pandas
pandas.options.mode.chained_assignment = None # default = "warn"
import networkx

# Custom
import partner.utility as putility
import partner.regression as pro_reg
import partner.decomposition as decomp
import stragglers.mbpdb_stratification as bpd_strat # problem when executing uk_biobank not as sub-directory...

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
    paths["mbpdb_regression"] = os.path.join(path_dock, "mbpdb_regression")

    # Remove previous files to avoid version or batch confusion.
    if restore:
        putility.remove_directory(path=paths["mbpdb_regression"])
    # Initialize directories.
    putility.create_directories(
        path=paths["mbpdb_regression"]
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
        path_dock, "mbpdb_organization",
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


def read_source_cohort_model_reference(
    path_dock=None,
    report=None,
):
    """
    Reads and organizes source information from file.

    arguments:
        path_dock (str): path to dock directory for source and product
            directories and files
        report (bool): whether to print reports

    raises:

    returns:
        (dict<object>): collection of Pandas data-frame tables with entry names
            (keys) derived from original names of files

    """

    # Define path to parent directory.
    path_directory_parent = os.path.join(
        path_dock, "parameters", "stragglers", "regression_cohorts_models",
    )
    # Read all files within parent directory and organize tables.
    pail = putility.read_all_pandas_tables_files_within_parent_directory(
        path_directory_parent=path_directory_parent,
        types_pandas_table_read={
            "execution": "int",
            "cohort": "string",
            "cohort_sort": "int",
            "dependence": "string",
            "dependence_sort": "int",
            "dependence_type": "string",
            "independence": "string",
            "model": "string",
            "model_sort": "int",
            "model_adjustment": "string",
            "model_context": "string",
            "model_note": "string",
        },
        report=report,
    )
    # Return information.
    return pail


##########
# Driver


def stratify_cohorts_call_run_regressions(
    table=None,
    table_cohorts_models=None,
    independences_summary=None,
    filter_execution=None,
    type=None,
    report=None,
):
    """
    Stratify cohorts for phenotype tables and call driver for regressions.

    Format of "table"...
    Pandas data frame with variables (features) across columns and
    their samples (cases, observations) across rows, with an explicit index

    Format of "table_cohorts_models"...
    columns: "execution", "cohort", "cohort_sort", "dependence",
    "dependence_sort", "dependence_type", "model", "model_sort", "independence",
    "model_note"

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        table_cohorts_models (object): Pandas data frame that specifies cohorts,
            dependent variables, and independent variables (models) for
            regression
        independences_summary (list<str>): names of independent variables for
            which to include information in the summary table, or "None" to
            include information for all original independent variables
        filter_execution (bool): whether to filter records in cohort-model table
            by logical binary "execution" variable
        type (str): type of regression analysis, either 'linear' or 'logistic'
        report (bool): whether to print reports

    raises:

    returns:
        (dict): information from regressions

    """

    # Stratify phenotypes in cohorts.
    records_cohorts = bpd_strat.stratify_phenotype_cohorts(
        table=table,
    )
    entries_cohorts = (
        putility.structure_from_records_to_entries(
            records=records_cohorts,
    ))
    # Call driver for regressions.
    pail = (
        pro_reg.drive_linear_logistic_regressions_cohorts_models(
            entries_cohorts=entries_cohorts,
            table_cohorts_models=table_cohorts_models,
            independences_summary=independences_summary,
            filter_execution=filter_execution,
            type=type,
            report=report,
    ))
    # Return information.
    return pail


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
    putility.print_terminal_partition(level=1)
    print(path_dock)
    print("version check: TCW, 17 April 2023")
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
    source_reference = read_source_cohort_model_reference(
        path_dock=path_dock,
        report=True,
    )

    # Collect information.
    pail_write = dict()
    pail_write["tables"] = dict()

    # Drive regressions.

    if True:
        pail_regression = stratify_cohorts_call_run_regressions(
            table=source["table_phenotypes"],
            table_cohorts_models=(
                source_reference["table_2023-04-24_tsh_polygenic_scores"]
            ),
            independences_summary=None, # "None" or list of variables
            filter_execution=True,
            type="linear",
            report=True,
        )
        pail_write["tables"]["table_2023-04-24_tsh_polygenic_scores"] = (
            pail_regression["table"]
        )

        pass

    # Write product information to file.
    pro_reg.write_product(
        pail_write=pail_write,
        path_directory=paths["mbpdb_regression"],
    )
    pass





#
