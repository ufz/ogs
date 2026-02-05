import itertools
import warnings
from pathlib import Path

from lxml import etree

"""
File that contains utility functions for the creation of the documentation and Feature Matrix generation.
"""


def getProjectFiles(path: Path) -> list[Path]:
    """Function to extract all file paths from OGS Test files"""
    path_list_xml = list(path.rglob("*.xml"))
    path_list = []

    for filepath in path_list_xml:
        with open(filepath, "r") as included_xml_file:
            if "OpenGeoSysProjectDiff" in included_xml_file.read():
                path_list.append(filepath)

    path_list.extend(list((path).rglob("*.prj")))
    return path_list


def getXMLFiles(path: Path):
    files = getProjectFiles(path=path)
    if not files:
        msg = "No project files found"
        raise ValueError(msg)

    # Parse the XML Files
    filenames = [file.name for file in files]
    xml_files = [xmlParser(file) for file in files]
    filenames = list(
        itertools.compress(filenames, [xml is not None for xml in xml_files])
    )
    files = list(itertools.compress(files, [xml is not None for xml in xml_files]))
    return [xml for xml in xml_files if xml is not None]


def xmlParser(filedir: Path) -> etree.ElementTree:
    """wrapper for etree.parse, as there are files present, that could not be parsed."""
    try:
        return etree.parse(filedir)
    except Exception:
        warnings.warn(f"Not able to parse: {filedir}.", stacklevel=2)
