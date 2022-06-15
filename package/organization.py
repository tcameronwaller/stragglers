
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
import bipolar_biobank.stratification as bpd_strat
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


def interpret_biological_genetic_sex_y(
    value_source=None,
):
    """
    Intepret representation of biological sex.

    The indication of genetic sex comes from a file in PLINK2 ".fam" format.
    https://www.cog-genomics.org/plink/2.0/formats#fam
    1: "male"
    2: "female"
    0: "unknown"

    arguments:
        value_source (str): raw value as character string

    raises:

    returns:
        (float): interpretation value (1: male, 0: female, NAN: missing, null)

    """

    # Interpret field code.
    if (
        (not pandas.isna(value_source)) and
        (len(str(value_source)) > 0)
    ):
        # The value is non-missing and hopefully interpretable.
        # Determine whether the value matches any strings.
        if (str(value_source) == "1"):
            # 1: "male"
            value_product = 1
        elif (str(value_source) == "2"):
            # 2: "female"
            value_product = 0
        elif (str(value_source) == "0"):
            # 0: "unknown"
            value_product = float("nan")
        else:
            # Ambiguous, uninterpretable, or missing information.
            value_product = float("nan")
    else:
        # Ambiguous, uninterpretable, or missing information.
        value_product = float("nan")
    # Return.
    return value_product



def determine_consensus_biological_sex_y(
    field_31=None,
    field_22001=None,
    field_22019=None,
):
    """
    Determine consensus biological sex (not social gender).

    Prioritize interpretation of the genetic sex variable, and use self report
    to fill in any missing values. Exclude records for persons whose genotypes
    suggested aneuploidy in the sex chromosomes.

    Use a logical binary encoding for presence of Y chromosome.
    0 : female (false, XX)
    1 : male (true, XY)

    Code is same as UK Biobank for data-fields "31" and "22001".
    0: "Female"
    1: "Male"

    arguments:
        field_31 (float): UK Biobank data-field 31, person's self-reported sex
        field_22001 (float): UK Biobank data-field 22001, genetic sex
        field_22019 (float): UK Biobank data-field 22019, sex-chromosome
            aneuploidy

    raises:

    returns:
        (float): interpretation value

    """

    # Interpret biological sex from genotype record.
    sex_y_genotype = interpret_genotype_biological_sex_y(
        value_source=sex_genetic_raw,
    )
    # Interpret biological sex from phenotype record.
    sex_y_phenotype = interpret_phenotype_biological_sex_y(
        value_source=,
    )
    # Comparison.
    # Determine whether there was aneuploidy in the sex chromosomes.
    if (
        (not pandas.isna(aneuploidy_sex)) and
        (aneuploidy_sex == 0)
    ):
        # Prioritize genetic sex.
        if (not pandas.isna(sex_genetic)):
            # Genetic sex variable has a valid value.
            # Prioritize this variable.
            value = sex_genetic
        elif (not pandas.isna(sex_self_report)):
            # Person has missing value for genetic sex.
            # Self-report sex variable has a valid value.
            # Resort to self-report sex in absence of genetic sex.
            value = sex_self_report
        else:
            # Missing information.
            value = float("nan")
    else:
        # Aneuploidy in the sex chromosomes complicates definition of biological
        # sex.
        value = float("nan")
    # Return information.
    return value




def interpret_biological_sex_text(
    sex_y=None,
):
    """
    Intepret textual representation of biological sex from logical binary
    representation of presence of Y chromosome.

    arguments:
        sex_y (float): logical binary representation of presence of Y chromosome

    raises:

    returns:
        (str): textual representation of biological sex

    """

    # Interpret field code.
    if (
        (not pandas.isna(sex_y)) and
        ((sex_y == 0) or (sex_y == 1))
    ):
        # The value is non-missing and hopefully interpretable.
        if (sex_y == 0):
            # sex_y: 0, "female"
            sex_text = "female"
        elif (sex_y == 1):
            # sex_y: 1, "male"
            sex_text = "male"
        else:
            # Ambiguous, uninterpretable, or missing information.
            sex_text = ""
    else:
        # Ambiguous, uninterpretable, or missing information.
        sex_text = ""
    # Return.
    return sex_text


def determine_biological_sex_variables(
    table=None,
    report=None,
):
    """
    Organizes table of information about phenotypes.

    arguments:
        table (object): Pandas data frame of information about phenotypes
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of information about phenotypes

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Interpret diagnosis type of bipolar disorder.
    table["sex_y"] = table.apply(
        lambda row:
            interpret_biological_genetic_sex_y(
                sex_genetic_raw=row["sex_genetic_raw"],
            ),
        axis="columns", # apply function to each row
    )
    table["sex_text"] = table.apply(
        lambda row:
            interpret_biological_sex_text(
                sex_y=row["sex_y"],
            ),
        axis="columns", # apply function to each row
    )

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        name_function = (
            "determine_biological_sex_variables()"
        )
        print(name_function)
        utility.print_terminal_partition(level=3)
        # Stratify tables.
        table_report_female = table.loc[
            (
                (table["sex_text"] == "female")
            ), :
        ]
        table_report_male = table.loc[
            (
                (table["sex_text"] == "male")
            ), :
        ]
        # Count.
        count_female = table_report_female.shape[0]
        count_male = table_report_male.shape[0]
        utility.print_terminal_partition(level=5)
        print("count female: " + str(count_female))
        utility.print_terminal_partition(level=5)
        print("count male: " + str(count_male))
        pass
    # Return information.
    return table






def interpret_bipolar_disorder_diagnosis_control_case(
    value_source=None,
    strings_case=None,
    strings_control=None,
):
    """
    Inteprets control or case of diagnosis of bipolar disorder.

    arguments:
        value_source (str): raw value as character string
        strings_case (list<str>): character strings that indicate diagnosis of
            Bipolar Disorder
        strings_control (list<str>): character strings that indiciate control

    raises:

    returns:
        (float): interpretation value (1: true, 0: false, NAN: missing, null)

    """

    if (
        (not pandas.isna(value_source)) and
        (len(str(value_source)) > 0)
    ):
        # The value is non-missing and hopefully interpretable.
        # Determine whether the value matches any strings.
        if (str(value_source) in strings_case):
            # 1: "case"
            value_product = 1
        elif (str(value_source) in strings_control):
            # 0: "control"
            value_product = 0
        else:
            # Assume that any other records are controls.
            value_product = 0
    else:
        # Some records for controls have empty strings.
        value_product = 0
    # Return.
    return value_product


def determine_logical_binary_indicator_variables_bipolar_disorder(
    table=None,
    report=None,
):
    """
    Organizes table of information about phenotypes.

    arguments:
        table (object): Pandas data frame of information about phenotypes
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of information about phenotypes

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Interpret diagnosis type of bipolar disorder.
    table["bipolar_disorder_control_case"] = table.apply(
        lambda row:
            interpret_bipolar_disorder_diagnosis_control_case(
                value_source=row["scid_dx"],
                strings_case=[
                    "BIPOLAR_I", "BIPOLAR_II",
                    "SCHIZOAFFECTIVE_BIPOLAR_TYPE", "OTHER",
                ],
                strings_control=["Not Eligible", "",],
            ),
        axis="columns", # apply function to each row
    )
    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        name_function = (
            "determine_logical_binary_indicator_variables_bipolar_disorder()"
        )
        print(name_function)
        utility.print_terminal_partition(level=3)
        # Stratify tables.
        table_report_control = table.loc[
            (
                (table["bipolar_disorder_control_case"] == 0)
            ), :
        ]
        table_report_case = table.loc[
            (
                (table["bipolar_disorder_control_case"] == 1)
            ), :
        ]
        # Count.
        count_control = table_report_control.shape[0]
        count_case = table_report_case.shape[0]
        utility.print_terminal_partition(level=5)
        print("count of Bipolar Disorder controls: " + str(count_control))
        utility.print_terminal_partition(level=5)
        print("count of Bipolar Disorder cases: " + str(count_case))
        pass
    # Return information.
    return table


def interpret_bipolar_disorder_type_diagnosis(
    value_source=None,
    match_1=None,
    match_0=None,
):
    """
    Inteprets whether the diagnosis type of bipolar disorder.

    arguments:
        value_source (str): raw value as character string
        match_1 (str): character string match for logical binary value one
        match_0 (str): character string match for logical binary value zero

    raises:

    returns:
        (float): interpretation value (1: true, 0: false, NAN: missing, null)

    """

    if (
        (not pandas.isna(value_source)) and
        (len(str(value_source)) > 0)
    ):
        # The value is non-missing and hopefully interpretable.
        # Determine whether the value matches any strings.
        if (str(value_source) == str(match_1)):
            # 1: "True"
            value_product = 1
        elif (str(value_source) == str(match_0)):
            # 0: "No"
            value_product = 0
        else:
            # Ambiguous, uninterpretable, or missing information.
            value_product = float("nan")
    else:
        # Ambiguous, uninterpretable, or missing information.
        value_product = float("nan")
    # Return.
    return value_product


def determine_logical_binary_indicator_variables_bipolar_disorder_type(
    table=None,
    report=None,
):
    """
    Organizes table of information about phenotypes.

    arguments:
        table (object): Pandas data frame of information about phenotypes
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of information about phenotypes

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Interpret diagnosis type of bipolar disorder.
    table["bipolar_disorder_type_1_2"] = table.apply(
        lambda row:
            interpret_bipolar_disorder_type_diagnosis(
                value_source=row["scid_dx"],
                match_1="BIPOLAR_I",
                match_0="BIPOLAR_II",
            ),
        axis="columns", # apply function to each row
    )
    table["bipolar_disorder_type_2_1"] = table.apply(
        lambda row:
            interpret_bipolar_disorder_type_diagnosis(
                value_source=row["scid_dx"],
                match_1="BIPOLAR_II",
                match_0="BIPOLAR_I",
            ),
        axis="columns", # apply function to each row
    )

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        name_function = (
            "determine_logical_binary_indicator_variables_bipolar_disorder_" +
            "type()"
        )
        print(name_function)
        utility.print_terminal_partition(level=3)
        # Stratify tables.
        table_report_type_1 = table.loc[
            (
                (table["bipolar_disorder_type_1_2"] == 1)
            ), :
        ]
        table_report_type_2 = table.loc[
            (
                (table["bipolar_disorder_type_2_1"] == 1)
            ), :
        ]
        # Count.
        count_type_1 = table_report_type_1.shape[0]
        count_type_2 = table_report_type_2.shape[0]
        utility.print_terminal_partition(level=5)
        print("count of Bipolar Disorder Type 1: " + str(count_type_1))
        utility.print_terminal_partition(level=5)
        print("count of Bipolar Disorder Type 2: " + str(count_type_2))
        pass
    # Return information.
    return table


def interpret_bipolar_disorder_rapid_cycling(
    value_source=None,
):
    """
    Intepret representation of rapid cycling as a phenotype of Bipolar Disorder.

    arguments:
        value_source (str): raw value as character string

    raises:

    returns:
        (float): interpretation value (1: True, 0: False, NAN: missing, null)

    """

    # Interpret field code.
    if (
        (not pandas.isna(value_source)) and
        (len(str(value_source)) > 0)
    ):
        # The value is non-missing and hopefully interpretable.
        # Determine whether the value matches any strings.
        if (str(value_source) == "1"):
            # 1: True
            value_product = 1
        elif (str(value_source) == "0"):
            # 0: "False"
            value_product = 0
        else:
            # Ambiguous, uninterpretable, or missing information.
            value_product = float("nan")
    else:
        # Ambiguous, uninterpretable, or missing information.
        value_product = float("nan")
    # Return.
    return value_product


def determine_logical_binary_indicator_variables_rapid_cycling(
    table=None,
    report=None,
):
    """
    Organizes table of information about phenotypes.

    arguments:
        table (object): Pandas data frame of information about phenotypes
        report (bool): whether to print reports

    raises:

    returns:
        (object): Pandas data frame of information about phenotypes

    """

    # Copy information in table.
    table = table.copy(deep=True)
    # Interpret diagnosis type of bipolar disorder.
    table["rapid_cycling"] = table.apply(
        lambda row:
            interpret_bipolar_disorder_rapid_cycling(
                value_source=row["rc"],
            ),
        axis="columns", # apply function to each row
    )

    # Report.
    if report:
        utility.print_terminal_partition(level=2)
        print("report: ")
        name_function = (
            "determine_logical_binary_indicator_variables_rapid_cycling()"
        )
        print(name_function)
        utility.print_terminal_partition(level=3)
        # Stratify tables.
        table_type_1_rapid_1 = table.loc[
            (
                (table["bipolar_disorder_type_1_2"] == 1) &
                (table["rapid_cycling"] == 1)
            ), :
        ]
        table_type_1_rapid_0 = table.loc[
            (
                (table["bipolar_disorder_type_1_2"] == 1) &
                (table["rapid_cycling"] == 0)
            ), :
        ]
        table_type_2_rapid_1 = table.loc[
            (
                (table["bipolar_disorder_type_2_1"] == 1) &
                (table["rapid_cycling"] == 1)
            ), :
        ]
        table_type_2_rapid_0 = table.loc[
            (
                (table["bipolar_disorder_type_2_1"] == 1) &
                (table["rapid_cycling"] == 0)
            ), :
        ]
        # Count.
        count_type_1_rapid_1 = table_type_1_rapid_1.shape[0]
        count_type_1_rapid_0 = table_type_1_rapid_0.shape[0]
        count_type_2_rapid_1 = table_type_2_rapid_1.shape[0]
        count_type_2_rapid_0 = table_type_2_rapid_0.shape[0]
        print(
            "count of Bipolar Disorder Type 1, Rapid Cycling True: " +
            str(count_type_1_rapid_1)
        )
        print(
            "count of Bipolar Disorder Type 1, Rapid Cycling False: " +
            str(count_type_1_rapid_0)
        )
        utility.print_terminal_partition(level=5)
        print(
            "count of Bipolar Disorder Type 2, Rapid Cycling True: " +
            str(count_type_2_rapid_1)
        )
        print(
            "count of Bipolar Disorder Type 2, Rapid Cycling False: " +
            str(count_type_2_rapid_0)
        )
        pass
    # Return information.
    return table


##########
# Write


def write_product_bipolar_organization(
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
    write_product_bipolar_organization(
        pail_write=pail_write["bipolar_organization"],
        path_directory=paths["bipolar_organization"],
    )
    pass




###############################################################################
# Procedure

# TODO: TCW; 14 June 2022
# TODO: interpret genetic sex and compare to gender
# TODO: interpret controls with neither Bipolar Disorder Type 1 nor 2
# TODO: 4. organize chronotype variables?


# TODO: use consensus of genetic sex and "gender" (prioritize genetic sex where available) for "sex_y"
# TODO: determine Bipolar Disorder controls and cases from genetic table


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
    # Determine variables for biological sex.
    table = determine_biological_sex_variables(
        table=source["table_phenotypes"],
        report=True,
    )
    # Determine logical binary indicator variables for Bipolar Disorder
    # diagnosis controls and cases.
    table = determine_logical_binary_indicator_variables_bipolar_disorder(
        table=table,
        report=True,
    )
    # Determine logical binary indicator variables for type of Bipolar Disorder
    # diagnosis.
    table = determine_logical_binary_indicator_variables_bipolar_disorder_type(
        table=table,
        report=True,
    )
    # Determine variables for rapid cycling as a phenotype of Bipolar Disorder.
    table = determine_logical_binary_indicator_variables_rapid_cycling(
        table=table,
        report=True,
    )

    # Stratify phenotype records in cohorts.
    bpd_strat.stratify_phenotype_cohorts(
        table=table,
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
        "identifier_genotype",
        "gwas1_sampleid",
        "gwas2_sampleid",
        "identifier_phenotype",
        "bib_id",
        "gender",
        "sex_genetic_raw",
        "sex_y",
        "sex_text",
        "control_case_raw",
        "bipolar_disorder_control_case",
        "pt_age",
        "BMI",
        "scid_dx",
        "bipolar_disorder_type_1_2",
        "bipolar_disorder_type_2_1",
        "rc",
        "rapid_cycling",
        "database",
        "steroid_globulin_female",
        "steroid_globulin_male",
        "testosterone_female",
        "testosterone_male",
    ]
    table = table.loc[
        :, table.columns.isin(columns_selection)
    ]
    table = table[[*columns_selection]]
    utility.print_terminal_partition(level=2)
    print("table after selection of columns")
    print(table)
    print("columns")
    print(table.columns.to_list())

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
