#!/usr/bin/python

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

# prevent broken pipe error
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

from print23 import print_
import os
import sys
import xml.etree.cElementTree as ET
import json

if len(sys.argv) != 3:
    sys.stderr.write("Usage:\n")
    sys.stderr.write("{0} DATADIR DOCAUXDIR\n".format(sys.argv[0]))
    sys.exit(1)

datadir = sys.argv[1]
docauxdir = sys.argv[2]

datadir = os.path.abspath(datadir)
docauxdir = os.path.abspath(docauxdir)
outdir = os.path.join(docauxdir, "dox", "CTestProjectFiles")

# Expansions or shortcuts for the documenation can be added here in the
# following format:
# "prj__processes__process": "process",
# See the expansion table in the append-xml-tags.py too.
tag_path_expansion_table = {
        "material__porous_medium__porous_medium" : "material__porous_medium",
        }

indent = "&nbsp;" * 2


def format_if_documented(is_doc, fmt, fullpagename, tag_attr, *args):
    if is_doc:
        tag_attr_formatted = r'\ref {0} "{1}"'.format(fullpagename, tag_attr)
    else:
        tag_attr_formatted = r'<span style="color: red;" title="undocumented: {0}">{1}</span>'.format(
            fullpagename, tag_attr)

    return fmt.format(tag_attr_formatted, tag_attr, *args)


def format_if_documented_nowarn(is_doc, fmt, fullpagename, tag_attr, *args):
    if is_doc:
        tag_attr_formatted = r'\ref {0} "{1}"'.format(fullpagename, tag_attr)
    else:
        tag_attr_formatted = tag_attr

    return fmt.format(tag_attr_formatted, tag_attr, *args)


def get_tagpath(pagename, typetag, typetag_levels_up):
    tagpath = pagename.replace("__", ".")
    is_doc = (tagpath, True) in documented_tags_attrs
    if not is_doc and typetag:
        tagpath_parts = tagpath.split(".")
        tagpath_parts.insert(-typetag_levels_up, typetag)
        tagpath_new = ".".join(tagpath_parts)
        is_doc = (tagpath_new, True) in documented_tags_attrs
        if is_doc:
            tagpath = tagpath_new
            pagename = tagpath.replace(".", "__")

    return tagpath, pagename, is_doc


def print_tags(node, level, pagename, fh, typetag, typetag_levels_up,
               relfilepath):
    global map_tag_to_prj_files
    tag = node.tag
    rootpagename = pagename

    if level > 0:
        pagename += "__" + tag if pagename else tag

    if pagename in tag_path_expansion_table:
        tagpath_unexpanded = pagename.replace("__", ".")
        dict_of_set_add(map_tag_to_prj_files, tagpath_unexpanded, relfilepath)
        pagename = tag_path_expansion_table[pagename]

    tagpath, pagename, is_doc = get_tagpath(pagename, typetag,
                                            typetag_levels_up)
    if is_doc:
        dict_of_set_add(map_tag_to_prj_files, tagpath, relfilepath)

    fullpagename = "ogs_file_param"
    if pagename:
        fullpagename += "__" + pagename

    attrs = ""
    for attr, value in sorted(node.attrib.items()):
        attrpath = tagpath + "." + attr
        attrpagename = "ogs_file_attr__" + attrpath.replace(".", "__")
        a_is_doc = (attrpath, False) in documented_tags_attrs
        attrs += format_if_documented(a_is_doc, ' {0}="{2}"', attrpagename,
                                      attr, value)
        if a_is_doc:
            dict_of_set_add(map_attr_to_prj_files, attrpath, relfilepath)

    if len(node) > 0:
        # node having child tag

        # print opening tag
        fh.write(format_if_documented(is_doc, \
                indent*level + '\\<{0}{2}\\><br>\n', fullpagename, tag, attrs))

        if node.text and node.text.strip():
            fh.write(indent * (level + 1) + node.text.strip() + "<br>\n")

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
            print_tags(child, level + 1, pagename, fh, typetag_children,
                       typetag_children_levels_up, relfilepath)

        # print closing tag
        fh.write((indent * level) + r"\</" + tag + "\\><br>\n")
    else:
        # node having no children
        if node.text and node.text.strip():
            if tag != "type":
                fh.write(format_if_documented(is_doc, \
                        indent*level + '\\<{0}{2}\\>{3}\\</{1}\\><br>\n', fullpagename, tag, attrs, node.text.strip()))
            else:
                typepagename = rootpagename + "__" + node.text.strip()
                typetagpath, typepagename, type_is_doc = get_tagpath(
                    typepagename, None, 0)
                if type_is_doc:
                    dict_of_set_add(map_tag_to_prj_files, typetagpath,
                                    relfilepath)
                typepagename = "ogs_file_param__" + typepagename

                # If the content of a type tag is undocumented no red
                # "undocumented..." HTML code will be generated.
                type_text_formatted = format_if_documented_nowarn(type_is_doc, \
                        '{0}', typepagename, node.text.strip())

                fh.write(format_if_documented(is_doc, \
                        indent*level + '\\<{0}{2}\\>{3}\\</{1}\\><br>\n', fullpagename, tag, attrs, type_text_formatted))
        else:
            fh.write(format_if_documented(is_doc, \
                    indent*level + '\\<{0}{2} /\\><br>\n', fullpagename, tag, attrs))


# set of (str tagpath, bool is_tag)
documented_tags_attrs = set()

# read parameter cache (generated by normalize-param-cache.py)
with open(os.path.join(docauxdir, "documented-parameters-cache.txt")) as fh:
    for line in fh:
        line = line.strip().split("@@@")
        if line[0] == "OK":
            tagpath = line[3]

            # assume parent tags are also documented
            tagpath_parts = tagpath.split(".")[:-1]
            while tagpath_parts:
                documented_tags_attrs.add((".".join(tagpath_parts), True))
                tagpath_parts.pop()

            method = line[-1]
            is_tag = method.find("Attribute") == -1
            documented_tags_attrs.add((tagpath, is_tag))
            # print_(tagpath, is_tag)


def has_prj_file_in_subdirs(reldirpath):
    for dn in dirs_with_prj_files:
        if dn.startswith(reldirpath):
            return True
    return False


dirs_with_prj_files = set()

# maps tags/attributes to the set of prj files they appear in
map_tag_to_prj_files = dict()
map_attr_to_prj_files = dict()


def dict_of_set_add(dos, key, value):
    if key not in dos: dos[key] = set()
    dos[key].add(value)


for (dirpath, dirnames, filenames) in os.walk(datadir, topdown=False):
    reldirpath = os.path.relpath(dirpath, datadir)
    outdirpath = os.path.join(outdir, reldirpath)
    print_(">", reldirpath)

    subpages = []

    for fn in filenames:
        filepath = os.path.join(dirpath, fn)
        relfilepath = os.path.relpath(filepath, datadir)
        pagename = "ogs_ctest_prj__" + relfilepath.replace("/", "__").replace(
            ".", "__")

        if fn.endswith(".prj"):
            outdoxfile = os.path.join(outdirpath, fn + ".dox")
            dirs_with_prj_files.add(reldirpath)

            subpages.append(pagename)

            if not os.path.exists(outdirpath):
                os.makedirs(outdirpath)

            with open(outdoxfile, "w") as fh:
                fh.write(r"""/*! \page %s %s

\parblock
<tt>
""" % (pagename, fn))

                try:
                    xmlroot = ET.parse(filepath).getroot()
                except ET.ParseError as err:
                    print("ParseError occured in file :", filepath)
                    print(err)
                    raise
                print_tags(xmlroot, 0, "prj", fh, None, 0, relfilepath)

                fh.write(r"""</tt>
\endparblock
*/
""")

    for fn in dirnames:
        filepath = os.path.join(dirpath, fn)
        relfilepath = os.path.relpath(filepath, datadir)
        if has_prj_file_in_subdirs(relfilepath):
            pagename = "ogs_ctest_prj__" + relfilepath.replace(
                "/", "__").replace(".", "__")
            subpages.append(pagename)

    if subpages:
        if not os.path.exists(outdirpath):
            os.makedirs(outdirpath)

        pagename = "ogs_ctest_prj__" + reldirpath.replace("/", "__").replace(
            ".", "__")
        pagetitle = os.path.split(reldirpath)[1]

        if pagetitle == ".":
            pagetitle = "OGS CTests&mdash;Project Files"

        with open(os.path.join(outdirpath, "index.dox"), "w") as fh:
            fh.write("""/*! \page {0} {1}

""".format(pagename, pagetitle))

            for sp in sorted(subpages):
                fh.write("- \\subpage {0}\n".format(sp))

            fh.write("""

*/""")

for k, v in map_tag_to_prj_files.items():
    map_tag_to_prj_files[k] = list(v)
    documented_tags_attrs.discard((k, True))

for k, v in map_attr_to_prj_files.items():
    map_attr_to_prj_files[k] = list(v)
    documented_tags_attrs.discard((k, False))

with open(os.path.join(docauxdir, "tested-parameters-cache.json"), "w") as fh:
    json.dump({ "tags": map_tag_to_prj_files, "attributes": map_attr_to_prj_files }, \
            fh, indent=2)

untested_tags = [tag for tag, is_tag in documented_tags_attrs if is_tag]
untested_attrs = [attr for attr, is_tag in documented_tags_attrs if not is_tag]

with open(os.path.join(docauxdir, "untested-parameters-cache.json"),
          "w") as fh:
    json.dump({ "tags": untested_tags, "attributes": untested_attrs }, \
            fh, indent=2)
