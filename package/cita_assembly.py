
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
import uk_biobank.organization as ukb_org
import promiscuity.utility as utility
import promiscuity.polygenic_score as pgs
#import promiscuity.plot as plot

###############################################################################
# Functionality

# TODO: 1. I need an "access" script for the CITA data set components
# TODO: 1.1. collect and organize paths to A. phenotypes B. polygenic scores



# TODO: Read in the Polygenic Scores (PGS) like I did for the Bipolar Biobank.
# TODO:
#    # Read and organize tables of polygenic scores.
#    tables_polygenic_scores = pgs.drive_read_organize_tables_polygenic_scores(
#        table_scores_parameters=source["table_scores_parameters"],
#        filter_inclusion=True,
#        report=True,
#    )




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
    path_table_parameter_scores = os.path.join(
        path_dock, "parameters", "bipolar_biobank",
        "polygenic_scores", "table_cita_biobank.tsv"
    )
    # Read information from file.
    table_parameter_scores = pgs.read_source_collection_polygenic_scores(
        path_table=path_table_parameter_scores,
        report=report,
    )
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
    # Collect and return information.
    return {
        "table_phenotypes": table_phenotypes,
        "table_parameter_scores": table_parameter_scores,
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


def convert_biochemical_concentration_units_moles_per_liter(
    table=None,
    factors_concentration=None,
):
    """
    Converts biochemical concentrations to units of moles per liter (mol/L).

    Molar Mass, Molecular Weight
    Species     ...     Molar Mass     ...     Reference
    estradiol           272.4 g/mol            PubChem
    testosterone        288.4 g/mol            PubChem
    SHBG
    albumin             69,367 g/mol           UniProt
    albumin             66,500 g/mol           Drug Bank
    - use 66.5 kDa molar mass for albumin (anticipate post-translational
    - - cleavage)

    Metric Prefixes
    (https://www.nist.gov/pml/weights-and-measures/metric-si-prefixes)
    Prefix     ...     Abbreviation     ...     Factor
    deci               d                        1E-1
    centi              c                        1E-2
    milli              m                        1E-3
    micro              u                        1E-6
    nano               n                        1E-9
    pico               p                        1E-12

    arguments:
        table (object): Pandas data frame table of phenotype variables
        factors_concentration (dict<float>): factors by which to multiply
            concentrations for convenience of storage and analysis

    raises:

    returns:
        (object): Pandas data frame table of phenotype variables

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Convert concentrations to units of moles per liter (mol/L).
    # Estradiol
    # original unit: pg / mL
    # final unit: pmol / L
    table["estradiol_pmol_l"] = table.apply(
        lambda row: float(
            (row["e2_"]) * (1E-12) * (1 / 272.4) * (1 / 1E-3)
            * factors_concentration["oestradiol"]
        ),
        axis="columns", # apply function to each row
    )
    # Testosterone
    # original unit: ng / dL
    # final unit: pmol / L
    table["testosterone_pmol_l"] = table.apply(
        lambda row: float(
            (row["testost"]) * (1E-9) * (1 / 288.4) * (1 / 1E-1)
            * factors_concentration["testosterone"]
        ),
        axis="columns", # apply function to each row
    )
    # Sex-Steroid Hormone Binding Globulin (SHBG)
    # original unit: nmol / L
    # final unit: nmol / L
    table["shbg_nmol_l"] = table.apply(
        lambda row: float(
            (row["shbg_"]) * (1E-9)
            * factors_concentration["steroid_globulin"]
        ),
        axis="columns", # apply function to each row
    )
    # Albumin
    # original unit: g / dL
    # final unit: umol / L
    table["albumin_umol_l"] = table.apply(
        lambda row: float(
            (row["albumin_"]) * (1 / 66500.0) * (1 / 1E-1)
            * factors_concentration["albumin"]
        ),
        axis="columns", # apply function to each row
    )
    # Return information.
    return table


def drive_calculation_estimate_bioavailable_free_hormones(
    factors_concentration=None,
    table=None,
    report=None,
):
    """
    Organizes calculation estimates of bioavailable and free concentrations
    of testosterone and oestradiol.

    arguments:
        factors_concentration (dict<float>): factors by which to multiply
            concentrations for the sake of float storage and analysis
        table (object): Pandas data frame table of phenotype variables
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame table of phenotype variables

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Define association constants for calculation of free hormones.
    # units: L/mol
    associations = dict()
    associations["shbg_test"] = 5.97E8 #
    associations["shbg_oest"] = 3.14E8 #
    associations["albumin_test"] = 4.06E4 #
    associations["albumin_oest"] = 4.21E4 #
    # Measurements.
    # Calculate estimation of free and bioavailable testosterone.
    table["testosterone_free_pmol_l"] = table.apply(
        lambda row:
            ukb_org.calculate_estimation_free_testosterone(
                testosterone=row["testosterone_pmol_l"],
                albumin=row["albumin_umol_l"],
                steroid_globulin=row["shbg_nmol_l"],
                factors_concentration=factors_concentration,
                associations=associations,
            ),
        axis="columns", # apply function to each row
    )
    table["testosterone_bioavailable_pmol_l"] = table.apply(
        lambda row:
            ukb_org.calculate_estimation_bioavailable_testosterone(
                testosterone_free=row["testosterone_free_pmol_l"],
                albumin=row["albumin_umol_l"],
                factors_concentration=factors_concentration,
                associations=associations,
            ),
        axis="columns", # apply function to each row
    )
    # Calculate estimation of free, bioavailable oestradiol.
    table["estradiol_free_pmol_l"] = table.apply(
        lambda row:
            ukb_org.calculate_estimation_free_oestradiol(
                oestradiol=row["estradiol_pmol_l"],
                testosterone_free=row["testosterone_free_pmol_l"],
                albumin=row["albumin_umol_l"],
                steroid_globulin=row["shbg_nmol_l"],
                factors_concentration=factors_concentration,
                associations=associations,
            ),
        axis="columns", # apply function to each row
    )
    table["estradiol_bioavailable_pmol_l"] = table.apply(
        lambda row:
            ukb_org.calculate_estimation_bioavailable_oestradiol(
                oestradiol=row["estradiol_pmol_l"],
                oestradiol_free=row["estradiol_free_pmol_l"],
                testosterone_free=row["testosterone_free_pmol_l"],
                steroid_globulin=row["shbg_nmol_l"],
                factors_concentration=factors_concentration,
                associations=associations,
            ),
        axis="columns", # apply function to each row
    )
    # Return information.
    return table


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
    table = convert_biochemical_concentration_units_moles_per_liter(
        table=table,
        factors_concentration=factors_concentration,
    )

    ##########
    # Calculate estimates of bioavailable and free hormones.
    table = drive_calculation_estimate_bioavailable_free_hormones(
        factors_concentration=factors_concentration,
        table=table,
        report=report,
    )

    # Return information.
    return table


def report_hormone_by_female_male(
    variable=None,
    table=None,
):
    """
    Reports counts and percentages of persons who were deficient in a hormone
    with stratification by sex, female menopause status, and age.

    arguments:
        column_hormone (str): name of column for hormone measurement
        threshold_deficiency (float): low threshold, concentrations below which
            qualify as deficiency
        table (object): Pandas data frame of phenotype variables across UK
            Biobank cohort

    raises:

    returns:

    """

    # Report.
    utility.print_terminal_partition(level=4)
    print("variable: " + str(variable))
    utility.print_terminal_partition(level=4)

    # Copy information in table.
    table = table.copy(deep=True)

    # Females.
    table_female = table.loc[
        (
            (table["gender"] == "Female")
        ), :
    ]
    table_male = table.loc[
        (
            (table["gender"] == "Male")
        ), :
    ]


    # Counts.
    count_female = int(table_female.shape[0])
    count_male = int(table_male.shape[0])
    # Quantitation.
    array_female = copy.deepcopy(table_female[variable].dropna().to_numpy())
    array_male = copy.deepcopy(table_male[variable].dropna().to_numpy())
    arrays = [array_female, array_male]
    for array in arrays:
        mean = numpy.nanmean(array)
        standard_error = scipy.stats.sem(array)
        interval_95 = (1.96 * standard_error)
        confidence_95_low = (mean - interval_95)
        confidence_95_high = (mean + interval_95)
        median = numpy.nanmedian(array)
        standard_deviation = numpy.nanstd(array)
        minimum = numpy.nanmin(array)
        maximum = numpy.nanmax(array)
        # Report.
        utility.print_terminal_partition(level=5)
        print("mean: " + str(round(mean, 3)))
        print("standard error: " + str(round(standard_error, 3)))
        print(
            "95% CI: " + str(
                str(round(confidence_95_low, 3)) + " ... " +
                str(round(confidence_95_high, 3))
            )
        )
    pass




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
    source = read_source(
        path_dock=path_dock,
        report=True,
    )
    # Read and organize tables of polygenic scores.
    tables_polygenic_scores = pgs.drive_read_organize_tables_polygenic_scores(
        table_parameter_scores=source["table_parameter_scores"],
        filter_inclusion=True,
        report=True,
    )




    if False:
        # This material belongs in "cita_organization" procedure.
        print(source["table_phenotypes"])

        # Calculate estimates of bioavailable and free estradiol and testosterone.
        table = drive_estimate_bioavailable_free_estradiol_testosterone(
            table=source["table_phenotypes"],
            report=True,
        )

        hormones = [
            "estradiol_pmol_l",
            "estradiol_bioavailable_pmol_l",
            "estradiol_free_pmol_l",
            "testosterone_pmol_l",
            "testosterone_bioavailable_pmol_l",
            "testosterone_free_pmol_l",
        ]
        for hormone in hormones:
            report_hormone_by_female_male(
                variable=hormone,
                table=table,
            )
    pass





#
