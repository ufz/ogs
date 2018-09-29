#!/usr/bin/python

# This script actually generates the QA page.
# For its usage see generate-project-file-doc-qa.sh

from print23 import print_
import sys
import re
import os.path
import json

github_src_url = "https://github.com/ufz/ogs/tree/master"


def debug(msg):
    sys.stderr.write(msg + "\n")


if len(sys.argv) != 3:
    print_("USAGE: {0} DOCAUXDIR SRCDIR".format(sys.argv[0]))
    sys.exit(1)

docauxdir = sys.argv[1]
if not os.path.isdir(docauxdir):
    print_("error: `{0}' is not a directory".format(docauxdir))
    sys.exit(1)

doxdir = os.path.join(docauxdir, "dox", "ProjectFile")
srcdir = sys.argv[2]

undocumented = []
unneeded_comments = []
wrong_input = []
no_doc_page = []
unneeded_md_files = dict()
good_tagpaths = set()
wrong_status = False

with open(os.path.join(docauxdir, "documented-parameters-cache.txt")) as fh:
    for inline in fh:
        inline = inline.strip().split("@@@")
        status = inline[0]

        if status == "OK":
            tag_path_comment = inline[3]
            tag_name_comment = tag_path_comment.split(".")[-1]

            dirs = tag_path_comment.split(".")[:-1]
            p = os.path.join(doxdir, *dirs)
            if os.path.isfile(os.path.join(p, "t_" + tag_name_comment + ".dox")) \
            or os.path.isfile(os.path.join(p, tag_name_comment, "i_" + tag_name_comment + ".dox")) \
            or os.path.isfile(os.path.join(p, tag_name_comment, "c_" + tag_name_comment + ".dox")) :
                good_tagpaths.add((tag_path_comment, "param"))
            elif os.path.isfile(
                    os.path.join(p, "a_" + tag_name_comment + ".dox")):
                good_tagpaths.add((tag_path_comment, "attr"))
            else:
                no_doc_page.append((tag_path_comment, inline[1], inline[2]))

        elif status == "WRONGIN":
            wrong_input.append(inline[1:])
        elif status == "NODOC":
            method = inline[6]
            # ignored parameters need not be documented
            if not method.startswith("ignore"):
                undocumented.append(inline[1:])
        elif status == "UNNEEDED":
            unneeded_comments.append(inline[1:])
        elif status == "SPECIAL":
            debug("SPECIAL: " +
                  " ".join(inline[1:]))  # TODO implement proper handling
            # unneeded.append(inline[1:])
        else:
            debug("ERROR: unrecognized status {0}".format(status))
            wrong_status = True

# traverse dox file hierarchy
srcdocdir = os.path.join(srcdir, "Documentation", "ProjectFile")
for (dirpath, _, filenames) in os.walk(srcdocdir):
    reldirpath = dirpath[len(srcdocdir) + 1:]

    for f in filenames:
        if not f.endswith(".md"): continue
        filepath = os.path.join(reldirpath, f)
        tag_or_attr = "param"

        if f.startswith("i_") or f.startswith("c_"):
            tagpath = reldirpath
        elif f.startswith("t_"):
            tagpath = os.path.join(reldirpath, f[2:-len(".md")])
        elif f.startswith("a_"):
            tagpath = os.path.join(reldirpath, f[2:-len(".md")])
            tag_or_attr = "attr"
        else:
            debug("ERROR: Found md file with unrecognized name: {0}"
                  .format(filepath))
            continue

        tagpath = tagpath.replace(os.sep, ".")

        if (tagpath, tag_or_attr) not in good_tagpaths:
            unneeded_md_files[(tagpath, tag_or_attr)] = filepath

qa_status_succeeded = True

# remove false positives from unneeded_md_files
if unneeded_md_files:
    qa_status_succeeded = False
    for tagpath, _ in good_tagpaths:
        tagpath = tagpath.split(".")
        while tagpath:
            tagpath.pop()
            parenttagpath = ".".join(tagpath)
            if (parenttagpath, "param") in unneeded_md_files:
                del unneeded_md_files[(parenttagpath, "param")]
                if not unneeded_md_files: break
        if not unneeded_md_files: break

if undocumented:
    qa_status_succeeded = False
    print_()
    print_("# Undocumented parameters")
    print_("| File | Line | Parameter | Type | Method | Link |")
    print_("| ---- | ---: | --------- | ---- | ------ | ---- |")
    for u in sorted(undocumented):
        u2 = list(u)
        u2.append(github_src_url)
        print_(("| {0} | {1} | {3} | <tt>{4}</tt> | <tt>{5}</tt> " +
                "| [&rarr; ufz/ogs/master]({6}/{0}#L{1})").format(*u2))

if unneeded_comments:
    qa_status_succeeded = False
    print_()
    print_("# Comments not documenting anything")
    print_("| File | Line | Comment | Link |")
    print_("| ---- | ---: | ------- | ---- |")
    for u in sorted(unneeded_comments):
        u2 = list(u)
        u2.append(github_src_url)
        u2[2] = re.sub(r'([\\@&$#<>%".|])', r"\\\1", u2[2])
        print_(("| {0} | {1} | {2} " +
                "| [&rarr; ufz/ogs/master]({3}/{0}#L{1}) |").format(*u2))

if wrong_input:
    qa_status_succeeded = False
    print_()
    print_("# Lines of input to that script that have not been recognized")
    print_("| File | Line | Content | Link |")
    print_("| ---- | ---: | ------- | ---- |")
    for w in sorted(wrong_input):
        w2 = list(w)
        w2.append(github_src_url)
        w2[2] = re.sub(r'([\\@&$#<>%".|])', r"\\\1", w2[2])
        print_(("| {0} | {1} | {2} " +
                "| [&rarr; ufz/ogs/master]({3}/{0}#L{1}) |").format(*w2))

if no_doc_page:
    qa_status_succeeded = False
    print_()
    print_("# No documentation page")
    print_("| Parameter | File | Line | Link |")
    print_("| --------- | ---- | ---: | ---- |")
    for n in sorted(no_doc_page):
        n2 = list(n)
        n2.append(github_src_url)
        print_(("| {0} | {1} | {2} " +
                "| [&rarr; ufz/ogs/master]({3}/{1}#L{2}) |").format(*n2))

if unneeded_md_files:
    qa_status_succeeded = False
    print_()
    print_("# Documentation pages that are not referenced in the source code")
    print_("| Page | *.md file | Link |")
    print_("| ---- | --------- | ---- |")
    for (tagpath, tag_or_attr), filepath in sorted(unneeded_md_files.items()):
        print_(
            (r'| \ref ogs_file_{0}__{1} | Documentation/ProjectFile/{2} ' +
             "| [&rarr; ufz/ogs/master]({3}/Documentation/ProjectFile/{2}#) |")
            .format(tag_or_attr,
                    tagpath.replace(".", "__"), filepath, github_src_url))

with open(os.path.join(docauxdir, "untested-parameters-cache.json")) as fh:
    untested_tags_attrs = json.load(fh)
    utags = [ ut for ut in untested_tags_attrs["tags"] \
            if ut != "gml" and not ut.startswith("gml.") ]
    if utags:
        qa_status_succeeded = False
        print_()
        print_("# Tags that do not occur in any CTest project file")
        for utag in sorted(utags):
            pagename = "ogs_file_param__" + utag.replace(".", "__")
            print_(r'- \ref {0} "{1}"'.format(pagename, utag))

    uattrs = [ ua for ua in untested_tags_attrs["attributes"] \
            if ua != "gml" and not ua.startswith("gml.") ]
    if uattrs:
        qa_status_succeeded = False
        print_()
        print_("# Attributes that do not occur in any CTest project file")
        for uattr in sorted(uattrs):
            pagename = "ogs_file_attr__" + uattr.replace(".", "__")
            print_(r'- \ref {0} "{1}"'.format(pagename, uattr))

# exit with error status if something was not documented/tested
if qa_status_succeeded:
    sys.exit(0)
sys.exit(1)
