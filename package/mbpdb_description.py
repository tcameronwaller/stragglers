"""
Organize information from the MBPDB project.

Author:

    T. Cameron Waller
    tcameronwaller@gmail.com
    Rochester, Minnesota 55904
    United States of America

License:

    This file is part of Stragglers
    (https://github.com/tcameronwaller/stragglers/).

    Stragglers supports analyses on data from multiple smaller sources.
    Copyright (C) 2022 Thomas Cameron Waller

    Stragglers is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the Free
    Software Foundation, either version 3 of the License, or (at your option)
    any later version.

    Stragglers is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
    more details.

    You should have received a copy of the GNU General Public License along
    with Stragglers. If not, see <http://www.gnu.org/licenses/>.
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
import stragglers.mbpdb_stratification as mb_strat
import promiscuity.utility as utility
import promiscuity.description as pdesc
#import promiscuity.plot as plot

###############################################################################
# Functionality

# TODO: TCW; 15 December 2022
# Introduce a new "attribution table" similar to the one in the UK Biobank project.
# Follow the pattern of counting various values of ordinal or categorical variables within a set of cohorts.


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
    paths["mbpdb_description"] = os.path.join(path_dock, "mbpdb_description")

    # Remove previous files to avoid version or batch confusion.
    if restore:
        utility.remove_directory(path=paths["mbpdb_description"])
    # Initialize directories.
    utility.create_directories(path=paths["mbpdb_description"])
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
# Table Attribution


def define_variables_table_attribution():
    """
    Defines values of nominal, categorical, or discrete variables for
    description in attribution table.

    arguments:

    raises:

    returns:
        (list<dict>): records with information about values of variables

    """

    # Collect records of information about each categorical varialbe and value.
    records = list()

    # General variables.

    # Variable: "sex_text"
    if False:

        record = dict()
        record["name"] = "sex_text_female"
        record["variable"] = "sex_text" # categorical or discrete variable
        record["value"] = "female" # categorical or discrete value of variable
        records.append(record)

        record = dict()
        record["name"] = "sex_text_male"
        record["variable"] = "sex_text" # categorical or discrete variable
        record["value"] = "male" # categorical or discrete value of variable
        records.append(record)

        pass

    # Variable: "bipolar_disorder_control_case"

    record = dict()
    record["name"] = "control"
    record["variable"] = "bipolar_disorder_control_case"
    record["value"] = 0
    records.append(record)

    record = dict()
    record["name"] = "case"
    record["variable"] = "bipolar_disorder_control_case"
    record["value"] = 1
    records.append(record)

    # Variable: "rapid_cycling"

    record = dict()
    record["name"] = "rapid_cycling_0"
    record["variable"] = "rapid_cycling"
    record["value"] = 0
    records.append(record)

    record = dict()
    record["name"] = "rapid_cycling_1"
    record["variable"] = "rapid_cycling"
    record["value"] = 1
    records.append(record)

    # Variable: "bipolar_disorder_type_1_2", "bipolar_disorder_type_2_1"

    record = dict()
    record["name"] = "bipolar_disorder_type_1"
    record["variable"] = "bipolar_disorder_type_1_2"
    record["value"] = 0
    records.append(record)

    record = dict()
    record["name"] = "bipolar_disorder_type_2"
    record["variable"] = "bipolar_disorder_type_2_1"
    record["value"] = 0
    records.append(record)


    # Return information
    return records


def organize_description_table_attribution(
    records_cohorts=None,
    report=None,
):
    """
    Organizes a description table for attribution of categorical or discrete
    variable values across cohorts.

    arguments:
        records_cohorts (list<dict>): records with information about cohorts
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of missingness of hormones in cohorts

    """

    # Define values of variables for attribution.
    records_attribution = define_variables_table_attribution()
    # Create a table from records of attribution.
    table_attribution = pdesc.drive_assemble_attribution_table(
        records_attribution=records_attribution,
        records_cohorts=records_cohorts,
        report=report,
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        print("organize_description_table_attribution()")
        utility.print_terminal_partition(level=3)
        pass
    # Return information.
    return table_attribution


##########
# Table Quantitation


def define_variables_table_quantitation():
    """
    Defines discrete or continuous variables on ordinal, interval, or ratio
    scales for description in quantitation table.

    arguments:

    raises:

    returns:
        (list<str>): names of variables

    """

    # Define variables.
    variables = [
        "age", "body", "genotype_pc_1",
    ]
    # Return information
    return variables


def organize_description_table_quantitation(
    records_cohorts=None,
    report=None,
):
    """
    Drives the assembly of a description table from records of quantitative
    descriptive statistics on variables of interest.

    These descriptive statistics are most appropriate for continuous variables
    on interval, or ratio scales, but they can also be informative for discrete
    variables on ordinal scales.

    arguments:
        records_cohorts (list<dict>): records with information about cohorts
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of missingness of hormones in cohorts

    """

    # Define variables.
    variables = define_variables_table_quantitation()
    # Create a table from records of quantitation.
    table_quantitation = pdesc.drive_assemble_quantitation_table(
        variables=variables,
        records_cohorts=records_cohorts,
        report=report,
    )

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        print("organize_description_table_quantitation()")
        utility.print_terminal_partition(level=3)
        pass
    # Return information.
    return table_quantitation



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

    # Stratify phenotype records in cohorts.
    records_cohorts = mb_strat.stratify_phenotype_cohorts(
        table=source["table_phenotypes"],
        report=True,
    )
    # Attribution table.
    table_attribution = organize_description_table_attribution(
        records_cohorts=records_cohorts,
        report=True,
    )
    #table_description = drive_collect_description_table_quantitation(
    #    variables=["age", "body", "genotype_pc_1",],
    #    records_cohorts=records_cohorts,
    #    report=True,
    #)

    # Quantitation table.
    table_quantitation = organize_description_table_quantitation(
        records_cohorts=records_cohorts,
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

    # Collect information.
    pail_write = dict()
    pail_write["table_attribution"] = table_attribution
    pail_write["table_quantitation"] = table_quantitation
    # Write product information to file.
    pdesc.write_product_tables(
        pail_write=pail_write,
        path_directory=paths["mbpdb_description"],
    )
    pass



#
