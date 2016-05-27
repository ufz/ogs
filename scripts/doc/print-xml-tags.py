#!/usr/bin/python

# prevent broken pipe error
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)

import os

import xml.etree.cElementTree as ET

import argparse

parser = argparse.ArgumentParser(description="Print XML tags")

parser.add_argument("ext",  help="Extension of files to consider")
parser.add_argument("path", help="Top level directory of traversal")

args = parser.parse_args()
rootdir = os.path.abspath(args.path)
extension = '.' + args.ext

# maps tags to the set of xml files they appear in
dict_tag_files = dict()

def dict_of_set_append(dict_, key, value):
    if key in dict_:
        dict_[key].add(value)
    else:
        dict_[key] = set((value,))


def print_tags(node, path, level, filepath):
    global dict_tag_files

    tag = node.tag
    if level>1: # skip root node
        tagpath = path + "." + tag
    else:
        tagpath = tag

    if level>0: # skip root node
        dict_of_set_append(dict_tag_files, "T | " + tagpath, filepath)
        for k in node.attrib:
            dict_of_set_append(dict_tag_files, "A | " + tagpath + "." + k, filepath)

    for child in node:
        print_tags(child, tagpath, level + 1, filepath)


for (dirpath, _, filenames) in os.walk(rootdir):
    for f in filenames:
        if not f.endswith(extension): continue

        filepath = os.path.join(dirpath, f)
        xmlroot = ET.parse(filepath).getroot()
        print_tags(xmlroot, "", 0, filepath[len(rootdir)+1:])

first = True
for (tag, files) in sorted(dict_tag_files.items()):
    if first:
        first = False
    else:
        print()

    print(tag)
    for f in sorted(files):
        print("   ", f)
