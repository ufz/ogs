#!/usr/bin/env python3

# This script augments the parameter documentation pages by information
# such as if they are required/optional, their data typ and in which
# end-to-end tests they are used.
# It uses the cache files generated by normalize-param-cache.py and by
# linked-xml-file.py

# prevent broken pipe error
from signal import SIG_DFL, SIGPIPE, signal

signal(SIGPIPE, SIG_DFL)

import json
import os
import sys

import pandas as pd
from print23 import print_

github_src_url = "https://gitlab.opengeosys.org/ogs/ogs/-/tree/master"
github_data_url = "https://gitlab.opengeosys.org/ogs/ogs/-/tree/master/Tests/Data"

if len(sys.argv) != 4:
    print_("Usage:")
    print_(f"{sys.argv[0]} EXT DATADIR DOCAUXDIR")
    sys.exit(1)

ext = sys.argv[1]
datadir = sys.argv[2]
docauxdir = sys.argv[3]

extension = "." + ext
datadir = os.path.abspath(datadir)
docauxdir = os.path.abspath(docauxdir)
docdir = os.path.join(docauxdir, "dox", "ProjectFile")

# used to expand documentation entry points to full xml tag paths
# that are used in the prj file.
# For example process on top-level could be expanded to processes.process with
#     "process":            "processes.process",
#
# See the expansion table in the linked-xml-file.py too.
tag_path_expansion_table = {
    "prj": "",
}


def write_parameter_type_info(fh, tagpath, tagpath_expanded, dict_tag_info):
    if tagpath:
        fh.write("\n\n# Additional info\n")
        if tagpath in dict_tag_info:
            for info in dict_tag_info[tagpath]:
                path = info[1]
                line = info[2]
                fh.write(f"\n## From {path} line {line}\n\n")

                method = info[6]

                try:
                    is_optional = info[7] == "True"
                except IndexError:
                    is_optional = False

                try:
                    is_defaulted = info[8] == "True"
                except IndexError:
                    is_defaulted = False

                try:
                    default_value = info[9]
                except IndexError:
                    default_value = ""

                if is_optional:
                    fh.write("- This is an optional parameter.\n")
                if is_defaulted:
                    fh.write(
                        f"- This parameter has a default value of <code>{default_value}</code>.\n"
                    )

                if method.endswith("List"):
                    fh.write("- This parameter can be given arbitrarily many times.\n")
                elif not is_optional:
                    fh.write("- This is a required parameter.\n")

                datatype = info[5]
                if datatype:
                    fh.write(f"- Data type: <code>{datatype}</code>\n")

                fh.write(f"- Expanded tag path: {tagpath_expanded}\n")

                fh.write(
                    "- Go to source code: [&rarr; ogs/ogs/master]({}/{}#L{})\n".format(
                        github_src_url, path, line
                    )
                )
        else:
            fh.write("\nNo additional info.\n")


def write_ctest_info(fh, tagpath, istag, tested_tags_attrs):
    fh.write("\n\n# Used in the following test data files\n\n")

    try:
        datafiles = tested_tags_attrs["tags" if istag else "attributes"][tagpath]
    except KeyError:
        fh.write("Used in no end-to-end test cases.\n")
        return

    for df in sorted(datafiles):
        pagename = "ogs_ctest_prj__" + df.replace("/", "__").replace(".", "__")
        fh.write(
            (
                "- \\[[&rarr; ogs/ogs/master]({1}/{0}) | "
                + '\\ref {2} "&rarr; doc"\\]&emsp;{0}\n'
            ).format(df, github_data_url, pagename)
        )


def write_ctest_media_info(fh, pcs_type, df_n_t_p_pcst):
    fh.write(
        """

# This process is commonly used together with the following media properties

Note: This list has been automatically extracted from OGS's benchmark tests (ctests).
Therefore it might not be exhaustive, but it should give users a good overview about
which properties they can/have to use with this process.
Probably most of the properties occurring in this list are mandatory.

The list might contain different property <tt>&lt;type&gt;</tt>s for some property
<tt>&lt;name&gt;</tt> to illustrate different possibilities the users have.

"""
    )

    df = df_n_t_p_pcst.set_index("pcs_type")

    try:
        df_this_pcs = df.xs(pcs_type)
    except KeyError:
        fh.write(
            f"""\
No media property information could be found in the ctests for the {pcs_type} process.
Either this process does not use <tt>&lt;media&gt;</tt> in the ctests or the extraction
scripts could not process the ctest prj files properly.

"""
        )
        return

    df_this_pcs_s = df_this_pcs.sort_values(["path", "name", "type"])

    old_path = None
    old_name = None

    for rec in df_this_pcs_s.itertuples(index=False):
        path = rec.path
        name = rec.name
        parts = path.split("/")

        if path != old_path:
            fh.write(
                f" - \\ref ogs_file_param__prj__{'__'.join(parts[1:])} \"{path}\"\n"
            )
            old_path = path
            old_name = None

        name_formatted = name
        if name != old_name:
            # highlight the first occurrence of each name in the list
            name_formatted = f"<b>{name}</b>"
            old_name = name

        type_ref = (
            ""
            if pd.isna(rec.type)
            else f"\\ref ogs_file_param__prj__{'__'.join(parts[1:])}__{rec.type}"
        )

        if pd.isna(rec.type):
            fh.write(f"   - <tt>&lt;name&gt;</tt> {name_formatted}\n")
        elif pd.isna(name):
            fh.write(f"   - <tt>&lt;type&gt;</tt> {type_ref}\n")
        else:
            fh.write(
                f'   - <div style="display: inline-block; min-width: 15em;"><tt>&lt;name&gt;</tt> {name_formatted}</div> <tt>&lt;type&gt;</tt> {type_ref}\n'
            )

    fh.write("\n\n")


def dict_of_list_append(dict_, key, value):
    if key in dict_:
        dict_[key].append(value)
    else:
        dict_[key] = [value]


# maps tags to additional parameter info obtained prior to this script
dict_tag_info = {}

with open(os.path.join(docauxdir, "tested-parameters-cache.json")) as fh:
    tested_tags_attrs = json.load(fh)

# read parameter cache (generated by normalize-param-cache.py)
with open(os.path.join(docauxdir, "documented-parameters-cache.txt")) as fh:
    for line in fh:
        line = line.strip().split("@@@")
        if line[0] == "OK":
            tagpath = line[3]
            dict_of_list_append(dict_tag_info, tagpath, line)

with open(os.path.join(docauxdir, "ctest-media-info.json")) as fh:
    df_n_t_p_pcst = pd.read_json(fh, orient="records")

# traverse dox file hierarchy
for dirpath, _, filenames in os.walk(docdir):
    reldirpath = dirpath[len(docdir) + 1 :]
    istag = True

    for f in filenames:
        if not f.endswith(".dox"):
            continue

        if f.startswith(("i_", "c_")):
            tagpath = reldirpath
        elif f.startswith("t_"):
            tagpath = os.path.join(reldirpath, f[2 : -len(".dox")])
            istag = True
        elif f.startswith("a_"):
            tagpath = os.path.join(reldirpath, f[2 : -len(".dox")])
            istag = False

        tagpath = tagpath.replace(os.sep, ".")

        path = os.path.join(dirpath, f)
        with open(path, "a") as fh:
            # TODO this can currently only expand the top level
            tagpathparts = tagpath.split(".")
            if tagpathparts[0] in tag_path_expansion_table:
                tagpathparts[0] = tag_path_expansion_table[tagpathparts[0]]
            else:
                tagpathparts[0] = "NONEXISTENT"
            tagpath_expanded = ".".join(tagpathparts).lstrip(".")

            # augment process documentation
            if (
                len(tagpathparts) > 0
                and tagpathparts[:-1]
                == [
                    "",  # "pcs" was replaced with empty string via the expansion table
                    "processes",
                    "process",
                ]
                and f.startswith("c_")
            ):
                pcs_type = tagpathparts[-1]
                write_ctest_media_info(fh, pcs_type, df_n_t_p_pcst)

            write_parameter_type_info(fh, tagpath, tagpath_expanded, dict_tag_info)

            if tagpath_expanded:
                write_ctest_info(fh, tagpath, istag, tested_tags_attrs)

            fh.write("\n*/\n")
