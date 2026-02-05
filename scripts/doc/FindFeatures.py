#!/usr/bin/env python


"""
Script and functions to find features of ogs and create a matrix which maps the respective features to
the Test project files.

This script is designed to extract information on features of the ogs test .prj files. This information
will also be written to the .json file, that is used to display the features on the html.

To add a feature: Add an element to the feature_dict in the get_feature_matrix function below.

To identify new features, there are several tools you can use. First the code_coverage attribute (or the
feature_matrix.getFilesWithLowestCodeCoverage function) of the Feature_matrix. It will show, what
percentage of code is already containing features. Be aware, that this can be misleading, as some of the
"unknown" features might already be part of a code block, that is connected to an existing feature.
You can also check the feature_matrix.getFilesWithLeastUsedFeatures Function to get project files, with
potential new features.
To have a nice way of showing, what lines are already connected to features use the
feature_matrix.getFileLines functions.
The functions: feature_matrix.getUnimportantFeatures1 and feature_matrix.getUnimportantFeatures2 are
useful to identify features that are tested sparsely.
"""


import argparse
import itertools
import json
import operator
from pathlib import Path
import re

from lxml import etree
import numpy as np
import pandas as pd

from FeatureMatrixClasses import FeatureMatrix, FeatureMatrixEntry
from utils import get_xml_files


def get_feature_dict(path: Path, xml_files: list[Path]) -> dict:
    # Create the feature dict
    feature_dict = {
        # Add Dummy features that are not "real" features but are needed to create the code coverage.
        "!Dummy: First Lines": lambda xml: checkFirstLines(xml),
        **{
            "!Dummy: "
            + case: lambda xml, case=case: check_tag_is_present(
                xml, "//" + case, line_type="open and close"
            )
            for case in [
                "OpenGeoSysProject",
                "media",
                "processes",
                "process",
                "phases",
                "medium",
                "properties",
                "time_loop",
                "nonlinear_solver",
                "parameters",
            ]
        },
        # Add Dummy feature for the Tags that should not appear on the website but used to calculate the code_coverage, so that it won't be considered as line without feature.
        "!Dummy: Comment": lambda xml: check_comment(xml),
        # Add Dummy feature for the Tags that should not appear on the website but used to calculate the code_coverage, so that it won't be considered as line without feature.
        "!Dummy: Media": lambda xml: check_tag_is_present(
            xml, "media", line_type="open and close"
        ),
        # Add Parameters. Checks whether there is a parameter type present. The possible types are extracted from the documentation.
        **{
            "Parameter: "
            + case: lambda xml, case=case: check_tag_text(
                xml, ".//parameters/parameter/type", case
            )
            for case in find_types_from_documentation(path, "parameters/parameter")
        },
        # Add Processes. Checks whether there is a process type present. The possible types are extracted from the documentation.
        **{
            "Process type: "
            + case: lambda xml, case=case: check_tag_text(
                xml, ".//processes/process/type", case
            )
            for case in find_types_from_documentation(path, "processes/process")
        },
        # Add Properties. Checks whether there is a process type present. The possible types are extracted from the documentation.
        **{
            "Property: "
            + case: lambda xml, case=case: check_tag_text(
                xml, ".//properties/property/type", case
            )
            for case in find_types_from_documentation(path, "../properties/property")
        },
        # Add Phases. Checks whether there are phase types present. The possible types are extracted from the documentation.
        **{
            "Phase: "
            + case: lambda xml, case=case: check_tag_text(
                xml, "./media/medium/phases/phase/type", case
            )
            for case in find_types_from_documentation(
                path, "./media/medium/phases/phase"
            )
        },
        # Add source terms for process variables. Checks whether there are source term types present. The possible types are extracted from the documentation.
        **{
            "Source_term: "
            + case: lambda xml, case=case: check_tag_text(
                xml,
                "./process_variables/process_variable/source_terms/source_term/type",
                case,
            )
            for case in find_types_from_documentation(
                path,
                "./process_variables/process_variable/source_terms/source_term",
            )
        },
        # Add boundary conditions for process variables. Checks whether there are boundary condition types present. The possible types are extracted from the documentation.
        **{
            "Boundary_condition: "
            + case: lambda xml, case=case: check_tag_text(
                xml,
                ".//process_variables/process_variable/boundary_conditions/boundary_condition/type",
                case,
            )
            for case in find_types_from_documentation(
                path,
                "./process_variables/process_variable/boundary_conditions/boundary_condition",
            )
        },
        # Add Convergence Criterion type for time_loops. The possible types are extracted from the documentation.
        **{
            "Convergence Criterion: "
            + case: lambda xml, case=case: check_tag_text(
                xml,
                "./time_loop/processes/process/convergence_criterion/type",
                case,
            )
            for case in find_types_from_documentation(
                path, "./time_loop/processes/process/convergence_criterion"
            )
        },
        # Add output type for time_stepping. The possible types are extracted from the documentation.
        **{
            "Time Stepping: "
            + case: lambda xml, case=case: check_tag_text(
                xml, "./time_loop/processes/process/time_stepping/type", case
            )
            for case in find_types_from_documentation(
                path, "./time_loop/processes/process/time_stepping"
            )
        },
        # Add output type for time_loops. The possible types are extracted from the documentation.
        **{
            "Output: "
            + case: lambda xml, case=case: check_tag_text(xml, "//output/type", case)
            for case in ["VTK", "XDMF"]
        },
        # Add multiple outputs.
        **{
            "Output: multiple": lambda xml: check_tag_children_count(
                xml, "//outputs", 1, ">"
            )
        },
        # Add order for process variables. Checks whether there are orders present. The possible entries are 1 or 2. Each of these will be given an entry in the feature matrix.
        **{
            "Order: "
            + case: lambda xml, case=case: check_tag_text(
                xml, "./process_variables/process_variable/order", case
            )
            for case in ["1", "2"]
        },
        # Add Component for process variables. Checks whether there are components present. The possible entries are 1, 2 or 3. Each of these will be given an entry in the feature matrix.
        **{
            "Component: "
            + case: lambda xml, case=case: check_tag_text(
                xml, "..//process_variables/process_variable/components", case
            )
            for case in ["1", "2", "3"]
        },
        # Check if there are linear solver tags with a child named: "eigen", "lis" or "petsc". For each of these a different function will be formulated.
        **{
            "Linear_solver: "
            + child: lambda xml, child_name=child: check_children_names(
                xml, ".//linear_solver", child_name
            )
            for child in ["eigen", "lis", "petsc"]
        },
        # Check whether there is a tag of the name given in the list. For each of the names a different Test-Function will be created.
        **{
            "Includes: "
            + tag_name: lambda xml, tag_name=tag_name: check_tag_is_present(
                xml, ".//" + tag_name
            )
            for tag_name in [
                "geometry",
                "curves",
                "local_coordinate_system",
                "python_script",
                "rasters",
                "meshes",
                "test_definition",
                "time_discretization",
                "global_process_coupling",
                "secondary_variable",
            ]
        },
        # Check for property names. For each possible name an feature dict entry will be created. The possible names are extracted by scanning all .prj files and look for said names.
        **{
            "Property name: "
            + name.replace("\n", ""): lambda xml, name=name: check_tag_text(
                xml, ".//properties/property/name", name
            )
            for name in np.unique(
                [
                    text
                    for xml in xml_files
                    for text in xml.xpath(".//properties/property/name/text()")
                ]
            )
        },  # All found entries for the names in all xml files.
        # Check for property names. For each possible name an feature dict entry will be created. The possible names are extracted by scanning all .prj files and look for said names.
        **{
            "Hdf n_files: "
            + name.replace("\n", ""): lambda xml, name=name: check_tag_text(
                xml, "//*/output/hdf/number_of_files", name
            )
            for name in np.unique(
                [
                    text
                    for xml in xml_files
                    for text in xml.xpath(".//*/output/hdf/number_of_files/text()")
                ]
            )
        },  # All found entries for the names in all xml files.
        # Check for Integration Order. For each possible name an feature dict entry will be created. The possible names are extracted by scanning all .prj files and look for said names.
        **{
            "Integration Order: "
            + name.replace("\n", ""): lambda xml, name=name: check_tag_text(
                xml, ".//processes/process/integration_order", name
            )
            for name in np.unique(
                [
                    text
                    for xml in xml_files
                    for text in xml.xpath(
                        ".//processes/process/integration_order/text()"
                    )
                ]
            )
        },  # All found entries for the names in all xml files.
        # Check whether there is a tag of the name of deactivated subdomain in the process_variables.
        "Includes: deactivated_subdomain ": lambda xml: check_tag_is_present(
            xml,
            ".//process_variables/process_variable/deactivated_subdomains",
        ),
        # Will check the given Tag and determine whether it is a tuple or a singular value.
        "Parameter Constant is Tuple:": lambda xml: check_tag_text_is_tuple(
            xml, '//parameters/parameter[type = "Constant"]/values'
        ),
        # Checks whether the Tag chemical System has an chemical solver attribute with the given values.
        **{
            "Chemical Solver:"
            + case: lambda xml, case=case: check_attributes(
                xml, "./chemical_system", "chemical_solver", case
            )
            for case in ["Phreeqc", "SelfContained"]
        },
        # Checks for damping of nonlinear solver.
        **{
            "Nonlinear solver:"
            + case: lambda xml, case=case: check_tag_is_present(
                xml, "//nonlinear_solvers/nonlinear_solver" + case
            )
            for case in ["damping", "damping_reduction", "recompute_jacobian"]
        },
        # Checks for type of nonlinear solver.
        **{
            "Nonlinear solver:"
            + case: lambda xml, case=case: check_tag_text(
                xml, "//nonlinear_solvers/nonlinear_solver/type", case
            )
            for case in ["Newton", "Picard", "PETScSNES"]
        },
        # Checks how many children meshes has.
        **{
            "Mesh Count: "
            + str(count): lambda xml, count=count: check_tag_children_count(
                xml, "./meshes", count
            )
            for count in np.unique(
                [
                    child_count
                    for xml in xml_files
                    for child_count in count_children(xml, "./meshes")
                ]
            )
        },
        # Check whether medium has an attribute id != 0
        "Medium: id not 0": lambda xml: check_attributes(
            xml,
            "./media/medium",
            "id",
            "0",
            operator="!=",
            line_type="open and close",
        ),
        # Check whether medium has an attribute id != 0
        "Medium: has phases": lambda xml: check_tag_is_present(
            xml, "./media/medium/phases", line_type="open and close"
        ),
        # Check whether mesh attribute axially symmetric is set to "true"
        "Mesh: axially_symmetric": lambda xml: check_attributes(
            xml, "mesh", "axially_symmetric", "true"
        ),
        # Check whether mesh attribute axially symmetric is set to "true"
        "Mesh: without <meshes>": lambda xml: check_tag_is_present(
            xml, "//OpenGeoSysProject/mesh"
        ),
        # Check if properties are present as children of the following "
        **{
            "Properties: child of "
            + case: lambda xml, case=case: check_tag_is_present(
                xml, "//" + case + "/properties"
            )
            for case in ["medium", "phase", "component"]
        },
        # Checks how many children meshes has.
        **{
            "Process Variables: "
            + str(count): lambda xml, count=count: check_tag_children_count(
                xml, "./process_variables", count
            )
            for count in np.unique(
                [
                    child_count
                    for xml in xml_files
                    for child_count in count_children(xml, "./process_variables")
                ]
            )
        },
        # Add submesh residuum output type for time_loops. The possible types are extracted from the documentation.
        **{
            "Submesh Residuum Output: "
            + case: lambda xml, case=case: check_tag_text(
                xml, "//submesh_residuum_output/type", case
            )
            for case in ["VTK", "XDMF"]
        },
        # check for Compensate Non Equilibrium Initial Residuum
        **{
            "Compensate Non Equilibrium Initial Residuum:"
            + case: lambda xml, case=case: check_tag_text(
                xml, "//compensate_non_equilibrium_initial_residuum", case
            )
            for case in ["true", "false"]
        },
    }
    return feature_dict


def get_feature_matrix(path: Path) -> FeatureMatrix:
    """
    First the project files are gathered. Afterwards the Feature Dictionaries are created. Eventually these dictionaries will be evaluated and
    the results will be put together in a matrix.
    The keys of the feature dict are the names of the features that will be displayed. The values should be functions that returns a
    FeatureMatrixEntry object. These functions start with "check*".
    """

    xml_files = get_xml_files(path, False)
    feature_dict = get_feature_dict(path, xml_files)

    return FeatureMatrix(feature_dict, xml_files)


########################################################################################
# Residual Functions                                               ####
########################################################################################


def create_json(mat: FeatureMatrix, filename: Path) -> str:
    """Creates JSON data from feature Matrix"""
    string = "[\n"

    for index, row in mat.has_feature.iterrows():
        string = (
            string
            + '  {\n    "file": "'
            + row.name.split("Tests/Data")[1].lstrip("/")
            + '",\n    "features": ['
        )
        string += ",".join(f'"{feature}"' for feature in mat.has_feature.columns[row])
        string += '],\n    "lines": ['
        intervals = extract_intervals(mat.lines, index)
        string += ",".join(f"{interval}" for interval in intervals if len(interval) > 0)

        string += '],\n    "feature_coverage": ' + f"{mat.code_coverage[index]}"
        string += "\n}"
        if index != mat.has_feature.index[-1]:
            string += ",\n"
        else:
            string += "\n"

    string += "]\n"
    with filename.open("w", encoding="utf-8") as text_file:
        text_file.write(string)
    return string


def check_attributes(
    xml: etree.ElementTree,
    xpath: str,
    attribute_name: str,
    case: str,
    line_type="range",
    operator: str = "==",
) -> FeatureMatrixEntry:
    """Checks whether the given Element has an attribute with the attribute_name given and will compare it using the operator given to the case attribute."""
    elements = xml.xpath(xpath)
    operator_func = translate_ops(operator)
    attributes = [
        operator_func(element.attrib[attribute_name], case)
        for element in elements
        if attribute_name in list(element.attrib.keys())
    ]
    return FeatureMatrixEntry(
        list(itertools.compress(elements, attributes)), line_type=line_type
    )


def check_children_names(
    xml: etree.ElementTree, xpath: str, case: str, line_type="range"
) -> FeatureMatrixEntry:
    """Will return a list of strings containing the children found in the xpath given along with the source lines of said child"""
    children = list(xml.xpath(xpath + "/*"))
    is_present = [case == child.tag for child in children]
    return FeatureMatrixEntry(list(itertools.compress(children, is_present)), line_type)


def check_comment(xml: etree.ElementTree) -> FeatureMatrixEntry:
    """Will find the lines, that contain comments in the xml files."""
    strings = etree.tostring(xml).split(b"\n")

    sourceline = xml.xpath("../OpenGeoSysProject")[
        0
    ].sourceline  # Number of lines before element starts

    from_lines = [
        line + sourceline - 1 for line in find_all_pattern_lines(strings, b"<!--")
    ]
    to_lines = [
        line + sourceline - 1 for line in find_all_pattern_lines(strings, b"-->")
    ]

    if len(from_lines) > 0:
        return FeatureMatrixEntry(
            [xml],
            lines=[
                pd.Interval(left=from_lines[i], right=to_lines[i], closed="both")
                for i in range(len(from_lines))
            ],
        )

    return FeatureMatrixEntry([])


def checkFirstLines(xml: etree.ElementTree) -> FeatureMatrixEntry:
    """Used as dummy feature to identify the lines before the <OpenGeoSysProject> Tag (prerequisites)."""
    xp = xml.xpath("../OpenGeoSysProject")

    if len(xp) > 0:
        if xp[0].sourceline - 1 < 1:
            return FeatureMatrixEntry(
                [xp[0]],
                lines=[
                    pd.Interval(
                        left=1,
                        right=1,
                        closed="both",
                    )
                ],
            )
        return FeatureMatrixEntry(
            [xp[0]],
            lines=[
                pd.Interval(
                    left=1,
                    right=xp[0].sourceline - 1,
                    closed="both",
                )
            ],
        )
    return FeatureMatrixEntry([])


def check_tag_children_count(
    xml: etree._ElementTree,
    xpath: str,
    count: int,
    operator: str = "==",
    line_type="range",
) -> FeatureMatrixEntry:
    """Will identify, how many children(Tags) this Tag has. Will compare the number using the operator to the count given. All Tags that match this comparison will be put in the FeatureMatrixEntry."""
    elements = xml.xpath(xpath)
    operator_func = translate_ops(operator)
    return FeatureMatrixEntry(
        list(
            itertools.compress(
                elements,
                [operator_func(len(el.getchildren()), count) for el in elements],
            )
        ),
        line_type,
    )


def check_tag_is_present(
    xml: etree.ElementTree, xpath: str, line_type="range"
) -> FeatureMatrixEntry:
    """Checks whether a tag is present in the parsed xml file under the given xpath. If True it will return the sourceline along with the boolean value"""
    return FeatureMatrixEntry(xml.xpath(xpath), line_type)


def check_tag_text(
    xml: etree.ElementTree,
    xpath: str,
    case: str,
    operator: str = "==",
    line_type="range",
) -> FeatureMatrixEntry:
    """Checks whether the text of the element of the xml file in the location of the given xpath is equal to the text given as case argument. If True will return the sourceline along with the boolean."""
    xpath_elements = xml.xpath(xpath)
    operator_func = translate_ops(operator)
    texts = [operator_func(x.text, case) for x in xpath_elements]

    if xpath.split("/")[-1] in ["name", "type"]:
        return FeatureMatrixEntry(
            [
                child.getparent()
                for child in list(itertools.compress(xpath_elements, texts))
            ],
            line_type,
        )
    return FeatureMatrixEntry(
        list(itertools.compress(xpath_elements, texts)), line_type
    )


def check_tag_text_is_tuple(
    xml: etree.ElementTree, xpath: str, line_type="range"
) -> FeatureMatrixEntry:
    """Checks whether the text of the element of the xml file in the location of the given xpath is a tuple. If True will return the sourceline along with the boolean."""
    elements = xml.xpath(xpath)
    is_tuple = [is_tuple_from_text(element.text) for element in elements]
    [element.text for element in elements]
    return FeatureMatrixEntry(list(itertools.compress(elements, is_tuple)), line_type)


def count_children(xml: etree._ElementTree, xpath: str) -> list[int]:
    # Counts the number of Children the tags gathered by the xpath have.
    return [len(element.getchildren()) for element in xml.xpath(xpath)]


def extract_intervals(
    lines: pd.DataFrame, index: object
) -> list[list[tuple[int, int]]]:
    intervals = []
    for f in range(len(lines.loc[index, :])):
        col_name = lines.columns[f]
        column_data = lines.loc[index, col_name]
        interval_list = [
            [int(interval.left), int(interval.right)] for interval in column_data
        ]
        intervals.append(interval_list)
    return intervals


def find_all_pattern_lines(xml_str: list[str], pattern: bytes) -> list[int]:
    return [
        j + 1
        for i in range(len(xml_str))
        if pattern in xml_str[i]
        for j in np.repeat(i, len(re.findall(pattern, xml_str[i])))
    ]


def find_types_from_documentation(path: Path, subdir: str) -> list[str]:
    # Will find all the possible types as defined in the documentation under the location "/Documentation/ProjectFile/prj/subdir/*"" and create a dictionary out of it
    base = path / "Documentation/ProjectFile/prj" / subdir
    return [
        m.group(1)
        for file in base.rglob("c_*.md")
        if (m := re.search(r"c_(.*)\.md", file.name))
    ]


def is_tuple_from_text(text: str) -> bool:
    """Evaluates Regex expression for given Text. If the text starts with a digit (possibly includes . or e) followed by a space and then another digit the
    function will return True else False."""
    return bool(
        re.compile(
            r"^\s*\d+(?:\.\d*)?(?:[eE][+-]?\d+)?\s+\d+(?:\.\d*)?(?:[eE][+-]?\d+)?\s*$"
        ).match(text)
    )


def translate_ops(op: str) -> operator:
    """ "Translates" operator strings into operator functions."""
    ops = {
        "!=": operator.__ne__,
        "in": operator.__contains__,
        ">": operator.gt,
        "<": operator.lt,
        ">=": operator.ge,
        "<=": operator.le,
        "==": operator.eq,
    }
    return ops[op]


def load_json_features(json_path: Path) -> list[dict]:
    """
    Loads a JSON file containing feature information for project files.

    The JSON file should have the following structure:
    [
      {
        "file": "path/to/file.prj",
        "features": ["feature1", "feature2", ...],
        "lines": [[[start1, end1], [start2, end2], ...], ...],
        "feature_coverage": 0.95
      },
      ...
    ]

    Args:
        json_path: Path to the JSON file

    Returns:
        List of dictionaries containing the parsed feature data
    """
    with Path.open(json_path, encoding="utf-8") as f:
        return json.load(f)


if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Generate feature matrix from directory"
    )
    parser.add_argument("--json", type=Path, help="Output JSON file name")
    parser.add_argument(
        "path",
        type=Path,
        metavar="PATH",
        nargs="?",
        default="./Tests/Data",
        help="Path to the directory containing feature files",
    )
    args = parser.parse_args()

    mat = get_feature_matrix(args.path)

    if args.json:
        create_json(mat, args.json)
        exit(0)
    create_json(mat, Path("./features.json"))
