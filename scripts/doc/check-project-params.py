#!/usr/bin/python

# This script actually generates the QA page.
# For its usage see generate-project-file-doc-qa.sh

import sys
import re
import os.path

github_src_url = "https://github.com/ufz/ogs/tree/master"

def debug(msg):
    sys.stderr.write(msg+"\n")

if len(sys.argv) != 2:
    print("USAGE: {} DOCAUXDIR".format(sys.argv[0]))
    sys.exit(1)

docauxdir = sys.argv[1]
if not os.path.isdir(docauxdir):
    print("error: `{}' is not a directory".format(docauxdir))
    sys.exit(1)

undocumented = []
unneeded_comments = []
wrong_input = []
no_doc_page = []

for inline in sys.stdin:
    inline = inline.strip().split("@@@")
    status = inline[0]

    if status == "OK":
        tag_path_comment = inline[3]
        tag_name_comment = tag_path_comment.split(".")[-1]

        dirs = tag_path_comment.split(".")[:-1]
        p = os.path.join(docauxdir, *dirs, )
        if     (not os.path.isfile(os.path.join(p,                   "t_" + tag_name_comment + ".dox"))) \
           and (not os.path.isfile(os.path.join(p,                   "a_" + tag_name_comment + ".dox"))) \
           and (not os.path.isfile(os.path.join(p, tag_name_comment, "i_" + tag_name_comment + ".dox"))) \
           and (not os.path.isfile(os.path.join(p, tag_name_comment, "c_" + tag_name_comment + ".dox"))) :
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
        debug("SPECIAL: " + " ".join(inline[1:])) # TODO implement proper handling
        # unneeded.append(inline[1:])
    else:
        debug("ERROR: unrecognized status {}".format(status))


if (undocumented):
    print()
    print("# Undocumented parameters")
    print("| File | Line | Parameter | Type | Method | Link |")
    print("| ---- | ---: | --------- | ---- | ------ | ---- |")
    for u in sorted(undocumented):
        print(("| {0} | {1} | {3} | <tt>{4}</tt> | <tt>{5}</tt> "
            + "| [&rarr; ufz/ogs/master]({6}/{0}#L{1})").format(*u, github_src_url))

if (unneeded_comments):
    print()
    print("# Comments not documenting anything")
    print("| File | Line | Comment | Link |")
    print("| ---- | ---: | ------- | ---- |")
    for u in sorted(unneeded_comments):
        u2 = list(u)
        u2[2] = re.sub(r'([\\@&$#<>%".|])', r"\\\1", u2[2])
        print(("| {0} | {1} | {2} "
            + "| [&rarr; ufz/ogs/master]({3}/{0}#L{1}) |").format(*u2, github_src_url))

if (wrong_input):
    print()
    print("# Lines of input to that script that have not been recognized")
    print("| File | Line | Content | Link |")
    print("| ---- | ---: | ------- | ---- |")
    for w in sorted(wrong_input):
        w2 = list(w)
        w2[2] = re.sub(r'([\\@&$#<>%".|])', r"\\\1", w2[2])
        print(("| {0} | {1} | {2} "
            + "| [&rarr; ufz/ogs/master]({3}/{0}#L{1}) |").format(*w2, github_src_url))

if (no_doc_page):
    print()
    print("# No documentation page")
    print("| Parameter | File | Line | Link |")
    print("| --------- | ---- | ---: | ---- |")
    for n in sorted(no_doc_page):
        print(("| {0} | {1} | {2} "
            + "| [&rarr; ufz/ogs/master]({3}/{1}#L{2}) |").format(*n, github_src_url))

# exit with error status if something was not documented.
if (not not undocumented) or (not not unneeded_comments) \
        or (not not wrong_input) or (not not no_doc_page):
            sys.exit(1)

sys.exit(0)
