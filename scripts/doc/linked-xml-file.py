#!/usr/bin/env python3

# This script traverses the DATADIR (where the ctest input data is) and creates
# out of every prj file found a file in a subdirectory of DOCAUXDIR.
#
# The created file is an input to Doxygen. The result of a subsequent Doxygen
# run will look the same as the original prj file (except for formatting), but
# each tag and attribute will be linked to the corresponding documentation page.
#
# A further output of this script are two JSON files, one telling which tag or
# attribute is tested in which prj file, and the other telling which tags and
# attributes are untested.

from collections import defaultdict
from pathlib import Path
from typing import TextIO
import xml.etree.ElementTree as ET
import sys
import os
import json
from signal import SIG_DFL, SIGPIPE, signal
import re

signal(SIGPIPE, SIG_DFL)


if len(sys.argv) != 3:
    sys.stderr.write("Usage:\n")
    sys.stderr.write(f"{sys.argv[0]} DATADIR DOCAUXDIR\n")
    sys.exit(1)

datadir = os.path.abspath(sys.argv[1])
docauxdir = os.path.abspath(sys.argv[2])
outdir = os.path.join(docauxdir, "dox", "CTestProjectFiles")
INDENT = "&nbsp;" * 2


def escape(s: str) -> str:
    """
    Escapes characters that could cause trouble with Doxygen.

    Also makes sure that newlines and indentation work.
    """

    # HTML escape sequence (&#47) does not work
    # <span> only works if an attribute (e.g. style) is present
    ESCAPE_SLASH = '<span style="">/</span>'

    # * " is treated specially sometimes, cf. https://www.doxygen.nl/manual/commands.html#cmdquot
    # * @ is special (JavaDoc commands)
    # * /* and */ start/end comments; we replace the / to avoid trouble with them
    s_escaped = s.replace('"', r"\"").replace("@", r"\@").replace("/", ESCAPE_SLASH)
    s_newline = s_escaped.replace("\n", "<br>\n")

    # finally, use &nbsp; for indentation
    return re.sub(
        r"^\s+",
        lambda m: "&nbsp;" * len(m.group(0)),
        s_newline,
        flags=re.MULTILINE,
    )


def format_if_documented(
    is_doc: bool, fmt: str, page_name: str, tag_attr: str, nowarn: bool, *args
) -> str:
    "Apply doxygen formatting to tag or attribute."

    tag_attr = escape(tag_attr)

    if is_doc:
        tag_attr_formatted = rf'\ref {page_name} "{tag_attr}"'
    elif nowarn:
        tag_attr_formatted = tag_attr
    else:
        tag_attr_formatted = (
            r'<span style="color: red;" title="undocumented: {}">{}</span>'.format(
                page_name, tag_attr
            )
        )
    return fmt.format(tag_attr_formatted, tag_attr, *args)


def replace_prefix(tag: str, typetag: str) -> str:
    "Replaces prefix to match the different docs structure for some tags."

    type_prefix = "" if typetag is None else typetag + "__"
    if (prefix := "prj__processes__process__constitutive_relation__") in tag:
        tag = tag.replace(
            prefix, "material.solid.constitutive_relation__" + type_prefix
        )

    for prefix_end in ["phases__phase__", "phases__", ""]:
        if (prefix := "prj__media__medium__" + prefix_end) in tag:
            if "property" not in tag:  # excludes parent tags
                return tag
            tag = tag.replace(prefix, "")

    return tag


def get_tagpath(
    page_name: str, typetag: str, typetag_levels_up: int
) -> tuple[str, str, bool]:
    tag_path = replace_prefix(page_name, typetag).replace("__", ".")
    is_doc = (tag_path, True) in documented_tags_attrs
    if not is_doc and typetag:
        tag_path_parts = tag_path.split(".")
        tag_path_parts.insert(-typetag_levels_up, typetag)
        tagpath_new = ".".join(tag_path_parts)
        is_doc = (tagpath_new, True) in documented_tags_attrs
        if is_doc:
            tag_path = tagpath_new
            page_name = tag_path.replace(".", "__")
    else:
        page_name = tag_path.replace(".", "__")

    return tag_path, page_name, is_doc


def remove_xml_header(xml_string: str) -> str:
    if xml_string.startswith("<?xml"):
        return "\n".join(xml_string.split("\n")[1:])
    return xml_string


def print_tags(
    node: ET.Element,
    level: int,
    page_name: str,
    filehandle: TextIO,
    typetag: str,
    typetag_levels_up: int,
    rel_filepath: str,
) -> None:
    tag = node.tag
    rootpagename = page_name

    # traverse down the included xml files
    if tag == "include":
        include_xml_path = (
            Path(datadir) / Path(rel_filepath).parent / node.attrib["file"]
        )
        with open(include_xml_path, "r") as included_xml_file:
            xml_string = remove_xml_header(included_xml_file.read())
        included_xml = ET.fromstring("<dummy>" + xml_string + "</dummy>")
        for included_child in included_xml:
            print_tags(
                included_child, level, page_name, filehandle, typetag,
                typetag_levels_up, rel_filepath
            )  # fmt: skip

    if level > 0:
        page_name += "__" + tag if page_name else tag

    tag_path, page_name, is_doc = get_tagpath(page_name, typetag, typetag_levels_up)

    if is_doc:
        map_tag_to_prj_files[tag_path].add(rel_filepath)

    fullpagename = "ogs_file_param"
    if page_name:
        fullpagename += "__" + page_name

    attrs = ""
    for attr, value in sorted(node.attrib.items()):
        attrpath = tag_path + "." + attr
        attrpagename = "ogs_file_attr__" + attrpath.replace(".", "__")
        a_is_doc = (attrpath, False) in documented_tags_attrs
        attrs += format_if_documented(
            a_is_doc,
            r" {0}=\"{2}\"",
            attrpagename,
            attr,
            True,
            escape(value),
        )
        if a_is_doc:
            map_attr_to_prj_files[attrpath].add(rel_filepath)

    if len(node) > 0:  # node having child tag
        fmt = INDENT * level + "\\<{0}{2}\\><br>\n"
        filehandle.write(
            format_if_documented(is_doc, fmt, fullpagename, tag, True, attrs)
        )

        if node.text and node.text.strip():
            filehandle.write(
                INDENT * (level + 1) + escape(node.text.strip()) + "<br>\n"
            )

        for child in node:
            if child.tag == "type" and child.text and child.text.strip():
                typetag_children = child.text.strip()
                typetag_children_levels_up = 1
                break
        else:
            typetag_children = typetag
            typetag_children_levels_up = typetag_levels_up + 1

        # Recursively process children tags.
        for child in node:
            print_tags(
                child,
                level + 1,
                page_name,
                filehandle,
                typetag_children,
                typetag_children_levels_up,
                rel_filepath,
            )

        # print closing tag
        filehandle.write((INDENT * level) + r"\</" + tag + "\\><br>\n")
    else:
        # node having no children
        if node.text and node.text.strip():
            if tag != "type":
                fmt = INDENT * level + "\\<{0}{2}\\>{3}\\</{1}\\><br>\n"
                filehandle.write(
                    format_if_documented(
                        is_doc, fmt, fullpagename, tag, True, attrs,
                        escape(node.text.strip()),
                    )  # fmt: skip
                )
            else:
                typepagename = rootpagename + "__" + node.text.strip()
                typetagpath, typepagename, type_is_doc = get_tagpath(
                    typepagename, None, 0
                )
                if type_is_doc:
                    map_tag_to_prj_files[typetagpath].add(rel_filepath)
                typepagename = "ogs_file_param__" + typepagename

                # If the content of a type tag is undocumented no red
                # "undocumented..." HTML code will be generated.
                type_text_formatted = format_if_documented(
                    type_is_doc,
                    "{0}",
                    typepagename,
                    escape(node.text.strip()),
                    False,
                )

                fmt = INDENT * level + "\\<{0}{2}\\>{3}\\</{1}\\><br>\n"
                filehandle.write(
                    format_if_documented(
                        is_doc, fmt, fullpagename, tag, True, attrs,
                        type_text_formatted,
                    )  # fmt: skip
                )
        else:
            fmt = INDENT * level + "\\<{0}{2} /\\><br>\n"
            filehandle.write(
                format_if_documented(is_doc, fmt, fullpagename, tag, True, attrs)
            )


# set of (str tagpath, bool is_tag)
documented_tags_attrs = set()

# read parameter cache (generated by normalize-param-cache.py)
with open(
    os.path.join(docauxdir, "documented-parameters-cache.txt"), encoding="UTF-8"
) as fh:
    for line in fh:
        line = line.strip().split("@@@")
        if line[0] == "OK":
            tagpath = line[3]

            # assume parent tags are also documented
            tagpath_parts = tagpath.split(".")[:-1]
            while tagpath_parts:
                documented_tags_attrs.add((".".join(tagpath_parts), True))
                tagpath_parts.pop()

            method = line[6]
            is_tag = method.find("Attribute") == -1
            documented_tags_attrs.add((tagpath, is_tag))


def has_prj_file_in_subdirs(rel_dirpath: str) -> bool:
    return any(dn.startswith(rel_dirpath) for dn in dirs_with_prj_files)


dirs_with_prj_files = set()

# maps tags/attributes to the set of prj files they appear in
map_tag_to_prj_files = defaultdict(set)
map_attr_to_prj_files = defaultdict(set)


for dirpath, dirnames, filenames in os.walk(datadir, topdown=False):
    reldirpath = os.path.relpath(dirpath, datadir)
    outdirpath = os.path.join(outdir, reldirpath)

    subpages = []

    for fn in filenames:
        filepath = os.path.join(dirpath, fn)
        relfilepath = os.path.relpath(filepath, datadir)
        pagename = "ogs_ctest_prj__" + relfilepath.replace("/", "__").replace(".", "__")

        if fn.endswith((".prj", ".xml")):
            if fn.endswith(".xml"):
                with open(filepath, "r") as included_xml_file:
                    if "OpenGeoSysProjectDiff" not in included_xml_file.read():
                        continue
            outdoxfile = os.path.join(outdirpath, fn + ".dox")
            dirs_with_prj_files.add(reldirpath)

            subpages.append(pagename)

            os.makedirs(outdirpath, exist_ok=True)

            # write linked prj file, cf. https://doxygen.opengeosys.org/d6/de3/ogs_ctest_prj__elliptic__circle_radius_1__circle_1e6_axi__prj
            with open(outdoxfile, "w", encoding="UTF-8") as fh:
                fh.write(
                    rf"""/*! \page {pagename} {fn}

\parblock
<tt>
"""
                )

                try:
                    xmlroot = ET.parse(filepath).getroot()
                except ET.ParseError as err:
                    print("ParseError occurred in file :", filepath)
                    print(err)
                    raise
                print_tags(xmlroot, 0, "prj", fh, None, 0, relfilepath)

                fh.write(
                    r"""</tt>
\endparblock
*/
"""
                )

    for fn in dirnames:
        filepath = os.path.join(dirpath, fn)
        relfilepath = os.path.relpath(filepath, datadir)
        if has_prj_file_in_subdirs(relfilepath):
            pagename = "ogs_ctest_prj__" + relfilepath.replace("/", "__").replace(
                ".", "__"
            )
            subpages.append(pagename)

    if subpages:
        os.makedirs(outdirpath, exist_ok=True)

        pagename = "ogs_ctest_prj__" + reldirpath.replace("/", "__").replace(".", "__")
        pagetitle = os.path.split(reldirpath)[1]

        if pagetitle == ".":
            pagetitle = "OGS CTests&mdash;Project Files"

        # Write CTest directory listing, cf. https://doxygen.opengeosys.org/de/d2a/ogs_ctest_prj__elliptic__circle_radius_1
        with open(os.path.join(outdirpath, "index.dox"), "w", encoding="UTF-8") as fh:
            fh.write(
                rf"""/*! \page {pagename} {pagetitle}

"""
            )

            for sp in sorted(subpages):
                fh.write(f"- \\subpage {sp}\n")

            fh.write(
                """

*/"""
            )

for k, v in map_tag_to_prj_files.items():
    map_tag_to_prj_files[k] = sorted(v)
    documented_tags_attrs.discard((k, True))

for k, v in map_attr_to_prj_files.items():
    map_attr_to_prj_files[k] = sorted(v)
    documented_tags_attrs.discard((k, False))

with open(
    os.path.join(docauxdir, "tested-parameters-cache.json"),
    "w",
    encoding="UTF-8",
) as fh:
    json.dump(
        {"tags": map_tag_to_prj_files, "attributes": map_attr_to_prj_files},
        fh,
        indent=2,
    )

untested_tags = sorted([tag for tag, is_tag in documented_tags_attrs if is_tag])
untested_attrs = sorted([attr for attr, is_tag in documented_tags_attrs if not is_tag])

# This goes to the QA page (https://doxygen.opengeosys.org/stable/d1/d49/project_file_doc_qa.html) :-)
with open(
    os.path.join(docauxdir, "untested-parameters-cache.json"),
    "w",
    encoding="UTF-8",
) as fh:
    json.dump({"tags": untested_tags, "attributes": untested_attrs}, fh, indent=2)
