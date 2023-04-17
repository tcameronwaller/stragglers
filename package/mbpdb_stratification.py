
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

# TODO: TCW; 15 December 2022
# I need new cohort stratifications for Rapid Cycling
# OR I need to implement at "attribution table" like I did in the UK Biobank project.


##########
# Stratification of cohorts


def stratify_phenotype_cohorts(
    table=None,
    report=None,
):
    """
    Organizes information and plots for sex hormones.

    arguments:
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort
        report (bool): whether to print reports

    raises:

    returns:
        (list<dict>): records with information about cohorts

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Collect records of information about each cohort and model.
    records = list()

    # Sex.

    record = dict()
    record["name"] = "female"
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female")
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "male"
    record["table"] = table.loc[
        (table["sex_text"] == "male"), :
    ]
    records.append(record)

    # Sex and genotype.

    record = dict()
    record["name"] = "genotype"
    record["table"] = table.loc[
        (
            (table["genotype_availability"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_genotype"
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["genotype_availability"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "male_genotype"
    record["table"] = table.loc[
        (
            (table["sex_text"] == "male") &
            (table["genotype_availability"] == 1)
        ), :
    ]
    records.append(record)

    # Sex and control.

    record = dict()
    record["name"] = "female_control"
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["bipolar_disorder_control_case"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "male_control"
    record["table"] = table.loc[
        (
            (table["sex_text"] == "male") &
            (table["bipolar_disorder_control_case"] == 0)
        ), :
    ]
    records.append(record)

    # Sex and Bipolar Disorder.

    record = dict()
    record["name"] = "female_bipolar_disorder"
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["bipolar_disorder_control_case"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "male_bipolar_disorder"
    record["table"] = table.loc[
        (
            (table["sex_text"] == "male") &
            (table["bipolar_disorder_control_case"] == 1)
        ), :
    ]
    records.append(record)

    # Sex and type of Bipolar Disorder.

    record = dict()
    record["name"] = "female_bipolar_disorder_1"
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["bipolar_disorder_type_1_2"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "male_bipolar_disorder_1"
    record["table"] = table.loc[
        (
            (table["sex_text"] == "male") &
            (table["bipolar_disorder_type_1_2"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "female_bipolar_disorder_2"
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["bipolar_disorder_type_2_1"] == 0)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "male_bipolar_disorder_2"
    record["table"] = table.loc[
        (
            (table["sex_text"] == "male") &
            (table["bipolar_disorder_type_2_1"] == 0)
        ), :
    ]
    records.append(record)

    # Sex and Rapid Cycling in Bipolar Disorder cases (any type).

    record = dict()
    record["name"] = "female_bipolar_disorder_rapid_cycle"
    record["table"] = table.loc[
        (
            (table["sex_text"] == "female") &
            (table["bipolar_disorder_control_case"] == 1) &
            (table["rapid_cycling"] == 1)
        ), :
    ]
    records.append(record)

    record = dict()
    record["name"] = "male_bipolar_disorder_rapid_cycle"
    record["table"] = table.loc[
        (
            (table["sex_text"] == "male") &
            (table["bipolar_disorder_control_case"] == 1) &
            (table["rapid_cycling"] == 1)
        ), :
    ]
    records.append(record)


    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        name_function = (
            "stratify_phenotype_cohorts()"
        )
        print(name_function)
        utility.print_terminal_partition(level=3)
        pass
    # Return information
    return records


###############################################################################
# Procedure
# Currently, this module is not executable.


#
