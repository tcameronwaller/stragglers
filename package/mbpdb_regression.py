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
import promiscuity.utility as utility
import promiscuity.regression as pro_reg
import promiscuity.decomposition as decomp
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
        utility.remove_directory(path=paths["mbpdb_regression"])
    # Initialize directories.
    utility.create_directories(
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
    pail = utility.read_all_pandas_tables_files_within_parent_directory(
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
        bpd_strat.organize_dictionary_entries_stratification_cohorts(
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


##########
# Write


def write_product_table(
    name=None,
    table=None,
    path_directory=None,
):
    """
    Writes product information to file.

    arguments:
        name (str): base name for file
        table (object): Pandas data-frame table to write to file
        path_directory (str): path to parent directory

    raises:

    returns:

    """

    # Reset index.
    table.reset_index(
        level=None,
        inplace=True,
        drop=True,
    )
    # Specify directories and files.
    path_table = os.path.join(
        path_directory, str(name + ".tsv")
    )
    # Write information to file.
    table.to_csv(
        path_or_buf=path_table,
        sep="\t",
        header=True,
        index=False,
    )
    pass


def write_product_tables(
    pail_write=None,
    path_directory=None,
):
    """
    Writes product information to file.

    arguments:
        pail_write (dict<object>): collection of information to write to file
        path_directory (str): path to parent directory

    raises:

    returns:

    """

    for name in pail_write.keys():
        write_product_table(
            name=name,
            table=pail_write[name],
            path_directory=path_directory,
        )
    pass


def write_product(
    pail_write=None,
    path_directory=None,
):
    """
    Writes product information to file.

    arguments:
        pail_write (dict<dict<object>>): collection of information to write to
            file
        path_directory (str): path to parent directory

    raises:

    returns:

    """

    # Export information.
    write_product_tables(
        pail_write=pail_write["tables"],
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
    print("version check: TCW, 18 August 2021")
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

    # Drive regressions.

    if True:
        pail_logistic_1 = stratify_cohorts_call_run_regressions(
            table=source["table_phenotypes"],
            table_cohorts_models=(
                source_reference["table_logistic_marginal_unadjust_bipolar_disorder_any_control_case"]
            ),
            independences_summary=None, # "None" or list of variables
            filter_execution=True,
            type="logistic",
            report=True,
        )
        pass
    if True:
        pail_logistic_2 = stratify_cohorts_call_run_regressions(
            table=source["table_phenotypes"],
            table_cohorts_models=(
                source_reference["table_logistic_joint_adjust_part_bipolar_disorder_any_control_case"]
            ),
            independences_summary=None, # "None" or list of variables
            filter_execution=True,
            type="logistic",
            report=True,
        )
        pass
    if True:
        pail_logistic_3 = stratify_cohorts_call_run_regressions(
            table=source["table_phenotypes"],
            table_cohorts_models=(
                source_reference["table_logistic_joint_adjust_full_bipolar_disorder_any_control_case"]
            ),
            independences_summary=None, # "None" or list of variables
            filter_execution=True,
            type="logistic",
            report=True,
        )
        pass



    if True:
        pail_logistic_4 = stratify_cohorts_call_run_regressions(
            table=source["table_phenotypes"],
            table_cohorts_models=(
                source_reference["table_logistic_marginal_unadjust_bipolar_disorder_any_rapid_cycling"]
            ),
            independences_summary=None, # "None" or list of variables
            filter_execution=True,
            type="logistic",
            report=True,
        )
        pass
    if True:
        pail_logistic_5 = stratify_cohorts_call_run_regressions(
            table=source["table_phenotypes"],
            table_cohorts_models=(
                source_reference["table_logistic_joint_adjust_part_bipolar_disorder_any_rapid_cycling"]
            ),
            independences_summary=None, # "None" or list of variables
            filter_execution=True,
            type="logistic",
            report=True,
        )
        pass
    if True:
        pail_logistic_6 = stratify_cohorts_call_run_regressions(
            table=source["table_phenotypes"],
            table_cohorts_models=(
                source_reference["table_logistic_joint_adjust_full_bipolar_disorder_any_rapid_cycling"]
            ),
            independences_summary=None, # "None" or list of variables
            filter_execution=True,
            type="logistic",
            report=True,
        )
        pass




    # Collect information.
    pail_write = dict()
    pail_write["tables"] = dict()

    # Logistic regression.

    pail_write["tables"]["table_logistic_marginal_unadjust_bipolar_disorder_any_control_case"] = (
        pail_logistic_1["table"]
    )
    pail_write["tables"]["table_logistic_joint_adjust_part_bipolar_disorder_any_control_case"] = (
        pail_logistic_2["table"]
    )
    pail_write["tables"]["table_logistic_joint_adjust_full_bipolar_disorder_any_control_case"] = (
        pail_logistic_3["table"]
    )

    pail_write["tables"]["table_logistic_marginal_unadjust_bipolar_disorder_any_rapid_cycling"] = (
        pail_logistic_4["table"]
    )
    pail_write["tables"]["table_logistic_joint_adjust_part_bipolar_disorder_any_rapid_cycling"] = (
        pail_logistic_5["table"]
    )
    pail_write["tables"]["table_logistic_joint_adjust_full_bipolar_disorder_any_rapid_cycling"] = (
        pail_logistic_6["table"]
    )

    # Linear regression.

    #pail_write["tables"]["table_bipolar_disorder_linear"] = (
    #    pail_linear_1["table"]
    #)
    # Write product information to file.
    write_product(
        pail_write=pail_write,
        path_directory=paths["mbpdb_regression"],
    )
    pass





#
