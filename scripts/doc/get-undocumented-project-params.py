#!/usr/bin/python

# expect input from get-project-params.sh

import sys
import subprocess

print_next = True

old_fn = None
undoc_lnos = []

def add_doc_stubs(fn, lnos):
    if not lnos: return

    print(fn, lnos)
    cmd = ["sed", "-i"]
    for lno in lnos:
        cmd.append("-e")
        cmd.append(str(lno) + r""" i \
//! \\ogs_file_param{todo_document_parameter} \\todo project_file_docu
//! \\ogs_file_param{todo_document_parameter} \\todo project_file_docu
""")
    cmd.append(fn)
    subprocess.run(cmd)
    del lnos[:]


for line in sys.stdin:
    fn, l, content = line.split(maxsplit=2)
    if fn != old_fn:
        add_doc_stubs(old_fn, undoc_lnos)
        old_fn = fn

    if content.startswith("//!"):
        print_next = False
    elif print_next:
        # print(line.rstrip())
        undoc_lnos.append(l)
    else:
        print_next = True

add_doc_stubs(old_fn, undoc_lnos)
