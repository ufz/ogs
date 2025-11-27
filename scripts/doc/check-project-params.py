#!/usr/bin/env python3

# This script actually generates the QA page.
# For its usage see generate-project-file-doc-qa.sh

import enum
import json
import os.path
import re
import sys
import traceback

from print23 import print_

github_src_url = "https://gitlab.opengeosys.org/ogs/ogs/-/tree/master"


class Status(enum.IntEnum):
    SUCCESS = 0
    WARN = 1
    FAIL = 2  # QA failure
    ERROR = 3  # exceptions, logic errors


def debug(msg):
    sys.stderr.write(msg + "\n")


def run():
    if len(sys.argv) != 3:
        print_(f"USAGE: {sys.argv[0]} DOCAUXDIR SRCDIR")
        sys.exit(1)

    docauxdir = sys.argv[1]
    if not os.path.isdir(docauxdir):
        print_(f"error: `{docauxdir}' is not a directory")
        sys.exit(1)

    doxdir = os.path.join(docauxdir, "dox", "ProjectFile")
    srcdir = sys.argv[2]

    undocumented = []
    unneeded_comments = []
    wrong_input = []
    no_doc_page = []
    unneeded_md_files = {}
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
                if (
                    os.path.isfile(os.path.join(p, "t_" + tag_name_comment + ".dox"))
                    or os.path.isfile(
                        os.path.join(
                            p, tag_name_comment, "i_" + tag_name_comment + ".dox"
                        )
                    )
                    or os.path.isfile(
                        os.path.join(
                            p, tag_name_comment, "c_" + tag_name_comment + ".dox"
                        )
                    )
                ):
                    good_tagpaths.add((tag_path_comment, "param"))
                elif os.path.isfile(os.path.join(p, "a_" + tag_name_comment + ".dox")):
                    good_tagpaths.add((tag_path_comment, "attr"))
                else:
                    no_doc_page.append((tag_path_comment, inline[1], inline[2]))

            elif status == "WRONGIN":
                wrong_input.append(inline[1:])
            elif status == "NODOC":
                method = inline[6]
                # ignored parameters need not be documented
                if not method.startswith("ignore"):
                    # after the method field there might be info about
                    # optional/default values, we ignore those
                    undocumented.append(inline[1:7])
            elif status == "UNNEEDED":
                unneeded_comments.append(inline[1:])
            elif status == "SPECIAL":
                debug(
                    "SPECIAL: " + " ".join(inline[1:])
                )  # TODO implement proper handling
                # unneeded.append(inline[1:])
            else:
                debug(f"ERROR: unrecognized status {status}")
                wrong_status = True

    # traverse dox file hierarchy
    srcdocdir = os.path.join(srcdir, "Documentation", "ProjectFile")
    for dirpath, _, filenames in os.walk(srcdocdir):
        reldirpath = dirpath[len(srcdocdir) + 1 :]

        for f in filenames:
            if not f.endswith(".md"):
                continue
            filepath = os.path.join(reldirpath, f)
            tag_or_attr = "param"

            if f.startswith(("i_", "c_")):
                tagpath = reldirpath
            elif f.startswith("t_"):
                tagpath = os.path.join(reldirpath, f[2 : -len(".md")])
            elif f.startswith("a_"):
                tagpath = os.path.join(reldirpath, f[2 : -len(".md")])
                tag_or_attr = "attr"
            else:
                debug(f"ERROR: Found md file with unrecognized name: {filepath}")
                continue

            tagpath = tagpath.replace(os.sep, ".")

            if (tagpath, tag_or_attr) not in good_tagpaths:
                unneeded_md_files[(tagpath, tag_or_attr)] = filepath

    qa_status = Status.SUCCESS

    # remove false positives from unneeded_md_files
    if unneeded_md_files:
        qa_status = False
        for tagpath, _ in good_tagpaths:
            tagpath = tagpath.split(".")
            while tagpath:
                tagpath.pop()
                parenttagpath = ".".join(tagpath)
                if (parenttagpath, "param") in unneeded_md_files:
                    del unneeded_md_files[(parenttagpath, "param")]
                    if not unneeded_md_files:
                        break
            if not unneeded_md_files:
                break

    if undocumented:
        qa_status = Status.FAIL
        print_()
        print_("# Undocumented parameters")
        print_("| File | Line | Parameter | Type | Method | Link |")
        print_("| ---- | ---: | --------- | ---- | ------ | ---- |")
        for u in sorted(undocumented):
            assert len(u) == 6
            print_(
                (
                    "| {0} | {1} | {3} | <tt>{4}</tt> | <tt>{5}</tt> "
                    + "| [&rarr; ogs/ogs/master]({6}/{0}#L{1})"
                ).format(*u, github_src_url)
            )
            print(
                "warning: undocumented parameter in {0}, line {1}: {3}".format(*u),
                file=sys.stderr,
            )

    if unneeded_comments:
        qa_status = Status.FAIL
        print_()
        print_("# Comments not documenting anything")
        print_("| File | Line | Comment | Link |")
        print_("| ---- | ---: | ------- | ---- |")
        for u in sorted(unneeded_comments):
            u2 = list(u)
            u2.append(github_src_url)
            u2[2] = re.sub(r'([\\@&$#<>%".|])', r"\\\1", u2[2])
            print_(
                (
                    "| {0} | {1} | {2} " + "| [&rarr; ogs/ogs/master]({3}/{0}#L{1}) |"
                ).format(*u2)
            )

    if wrong_input:
        qa_status = Status.FAIL
        print_()
        print_("# Lines of input to that script that have not been recognized")
        print_("| File | Line | Content | Link |")
        print_("| ---- | ---: | ------- | ---- |")
        for w in sorted(wrong_input):
            w2 = list(w)
            w2.append(github_src_url)
            w2[2] = re.sub(r'([\\@&$#<>%".|])', r"\\\1", w2[2])
            print_(
                (
                    "| {0} | {1} | {2} " + "| [&rarr; ogs/ogs/master]({3}/{0}#L{1}) |"
                ).format(*w2)
            )

    if no_doc_page:
        qa_status = Status.FAIL
        print_()
        print_("# No documentation page")
        print_("| Parameter | File | Line | Link |")
        print_("| --------- | ---- | ---: | ---- |")
        for n in sorted(no_doc_page):
            n2 = list(n)
            n2.append(github_src_url)
            print_(
                (
                    "| {0} | {1} | {2} " + "| [&rarr; ogs/ogs/master]({3}/{1}#L{2}) |"
                ).format(*n2)
            )

    if unneeded_md_files:
        qa_status = Status.FAIL
        print_()
        print_("# Documentation pages that are not referenced in the source code")
        print_("| Page | *.md file | Link |")
        print_("| ---- | --------- | ---- |")
        for (tagpath, tag_or_attr), filepath in sorted(unneeded_md_files.items()):
            print_(
                (
                    r"| \ref ogs_file_{0}__{1} | Documentation/ProjectFile/{2} "
                    + "| [&rarr; ogs/ogs/master]({3}/Documentation/ProjectFile/{2}#) |"
                ).format(
                    tag_or_attr, tagpath.replace(".", "__"), filepath, github_src_url
                )
            )

    with open(os.path.join(docauxdir, "untested-parameters-cache.json")) as fh:
        untested_tags_attrs = json.load(fh)
        utags = [
            ut
            for ut in untested_tags_attrs["tags"]
            if ut != "gml" and not ut.startswith("gml.")
        ]
        if utags:
            qa_status = max(qa_status, Status.WARN)
            print_()
            print_("# Tags that do not occur in any CTest project file")
            for utag in sorted(utags):
                pagename = "ogs_file_param__" + utag.replace(".", "__")
                print_(rf'- \ref {pagename} "{utag}"')

        uattrs = [
            ua
            for ua in untested_tags_attrs["attributes"]
            if ua != "gml" and not ua.startswith("gml.")
        ]
        if uattrs:
            qa_status = max(qa_status, Status.WARN)
            print_()
            print_("# Attributes that do not occur in any CTest project file")
            for uattr in sorted(uattrs):
                pagename = "ogs_file_attr__" + uattr.replace(".", "__")
                print_(rf'- \ref {pagename} "{uattr}"')

    return qa_status


if __name__ == "__main__":
    success = False
    try:
        qa_status = run()
    except:
        print(traceback.format_exc(), file=sys.stderr)
        sys.exit(Status.ERROR)

    sys.exit(qa_status)
