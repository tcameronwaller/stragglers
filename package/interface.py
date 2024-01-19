"""
Organizes parameters and manages execution of procedures within the 'stragglers'
package.

Execution of this subpackage, 'stragglers' occurs under the management of a
higher level package. Importation paths represent this hierarchy.

Author:

    T. Cameron Waller
    tcameronwaller@gmail.com
    Rochester, Minnesota 55904
    United States of America

License:

    This file is part of the 'stragglers' project
    (https://github.com/tcameronwaller/stragglers/).

    Project 'stragglers' supports processes and analyses on data from smaller
    studies and data sets.
    Copyright (C) 2022 Thomas Cameron Waller

    Project 'stragglers' is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    Project 'stragglers' is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
    more details.

    You should have received a copy of the GNU General Public License along
    with project 'stragglers'. If not, see <http://www.gnu.org/licenses/>.
"""
###############################################################################
# Notes

###############################################################################
# Installation and importation

# Standard.
import argparse
import textwrap

# Relevant.

# Custom.
import stragglers.scrap
import stragglers.mbpdb_prioritize_supplement
import stragglers.mbpdb_assembly
import stragglers.mbpdb_organization
import stragglers.mbpdb_description
import stragglers.mbpdb_regression
import stragglers.mcita_assembly
import stragglers.mcita_organization
#import stragglers.mcita_regression

#dir()
#importlib.reload()

###############################################################################
# Functionality


# Package's main subparser.


def define_subparser_main(subparsers=None):
    """
    Defines subparser and parameters.

    arguments:
        subparsers (object): reference to subparsers' container

    raises:

    returns:
        (object): reference to subparser

    """

    # Define description.
    description = stragglers.interface.define_description_main()
    # Define epilog.
    epilog = stragglers.interface.define_epilog_main()
    # Define parser.
    parser = subparsers.add_parser(
        name="stragglers",
        description=description,
        epilog=epilog,
        help="Help for 'stragglers' routine.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    # Define arguments.
    parser.add_argument(
        "-path_dock", "--path_dock", dest="path_dock", type=str, required=True,
        help=(
            "Path to dock directory for source and product " +
            "directories and files."
        )
    )
    parser.add_argument(
        "-scrap", "--scrap", dest="scrap",
        action="store_true",
        help=(
            "Scrap functionality for experimentation."
        )
    )
    parser.add_argument(
        "-mcita_assembly", "--mcita_assembly", dest="mcita_assembly",
        action="store_true",
        help=(
            "Assemble information."
        )
    )
    parser.add_argument(
        "-mcita_organization", "--mcita_organization",
        dest="mcita_organization", action="store_true",
        help=(
            "Organize information."
        )
    )
    parser.add_argument(
        "-mbpdb_prioritize_supplement", "--mbpdb_prioritize_supplement",
        dest="mbpdb_prioritize_supplement",
        action="store_true",
        help=(
            "Prioritize source of supplemental information."
        )
    )
    parser.add_argument(
        "-mbpdb_assembly", "--mbpdb_assembly", dest="mbpdb_assembly",
        action="store_true",
        help=(
            "Assemble information."
        )
    )
    parser.add_argument(
        "-mbpdb_organization", "--mbpdb_organization",
        dest="mbpdb_organization",
        action="store_true",
        help=(
            "Organize information beyond assembly procedure."
        )
    )
    parser.add_argument(
        "-mbpdb_description", "--mbpdb_description",
        dest="mbpdb_description",
        action="store_true",
        help=(
            "Describe information after organization."
        )
    )
    parser.add_argument(
        "-mbpdb_stratification", "--mbpdb_stratification",
        dest="mbpdb_stratification",
        action="store_true",
        help=(
            "Stratify cohorts and phenotypes for analyses."
        )
    )
    parser.add_argument(
        "-mbpdb_regression", "--mbpdb_regression",
        dest="mbpdb_regression",
        action="store_true",
        help=(
            "Regress phenotypes within cohorts."
        )
    )
    # Define behavior.
    parser.set_defaults(func=stragglers.interface.evaluate_parameters_main)
    # Return parser.
    return parser


def define_description_main():
    """
    Defines description for terminal interface.

    arguments:

    raises:

    returns:
        (str): description for terminal interface

    """

    description = textwrap.dedent("""\
        --------------------------------------------------
        --------------------------------------------------
        --------------------------------------------------

        Package's main procedure

        Do stuff.

        --------------------------------------------------
    """)
    return description


def define_epilog_main():
    """
    Defines epilog for terminal interface.

    arguments:

    raises:

    returns:
        (str): epilog for terminal interface

    """

    epilog = textwrap.dedent("""\

        --------------------------------------------------
        main routine

        --------------------------------------------------
        additional notes...


        --------------------------------------------------
        --------------------------------------------------
        --------------------------------------------------
    """)
    return epilog


def evaluate_parameters_main(arguments):
    """
    Evaluates parameters for model procedure.

    arguments:
        arguments (object): arguments from terminal

    raises:

    returns:

    """

    print("--------------------------------------------------")
    print("... call to main routine ...")
    # Execute procedure.
    if arguments.scrap:
        # Report status.
        print("... executing 'scrap' procedure ...")
        # Execute procedure.
        stragglers.scrap.execute_procedure(
            path_dock=arguments.path_dock
        )
    if arguments.mcita_assembly:
        # Report status.
        print("... executing 'mcita_assembly' procedure ...")
        # Execute procedure.
        stragglers.mcita_assembly.execute_procedure(
            path_dock=arguments.path_dock
        )
    if arguments.mcita_organization:
        # Report status.
        print("... executing 'mcita_organization' procedure ...")
        # Execute procedure.
        stragglers.mcita_organization.execute_procedure(
            path_dock=arguments.path_dock
        )
    if arguments.mbpdb_prioritize_supplement:
        # Report status.
        print("... executing 'mbpdb_prioritize_supplement' procedure ...")
        # Execute procedure.
        stragglers.mbpdb_prioritize_supplement.execute_procedure(
            path_dock=arguments.path_dock
        )
    if arguments.mbpdb_assembly:
        # Report status.
        print("... executing 'mbpdb_assembly' procedure ...")
        # Execute procedure.
        stragglers.mbpdb_assembly.execute_procedure(
            path_dock=arguments.path_dock
        )
    if arguments.mbpdb_organization:
        # Report status.
        print("... executing 'mbpdb_organization' procedure ...")
        # Execute procedure.
        stragglers.mbpdb_organization.execute_procedure(
            path_dock=arguments.path_dock
        )
    if arguments.mbpdb_description:
        # Report status.
        print("... executing 'mbpdb_description' procedure ...")
        # Execute procedure.
        stragglers.mbpdb_description.execute_procedure(
            path_dock=arguments.path_dock
        )
    if arguments.mbpdb_stratification:
        # Report status.
        print("... executing 'mbpdb_stratification' procedure ...")
        # Execute procedure.
        stragglers.mbpdb_stratification.execute_procedure(
            path_dock=arguments.path_dock
        )
    if arguments.mbpdb_regression:
        # Report status.
        print("... executing 'mbpdb_regression' procedure ...")
        # Execute procedure.
        stragglers.mbpdb_regression.execute_procedure(
            path_dock=arguments.path_dock
        )
    pass



###############################################################################
# Procedure



#
