#!/usr/bin/python

# This script takes the output of get-project-params.sh on stdin
# and transforms it into a tabular representation for further
# processing.

from print23 import print_
import sys
import re
import os.path

def debug(msg):
    sys.stderr.write(msg+"\n")

def write_out(*args):
    print_("@@@".join(str(a) for a in args))

# capture #2 is the parameter path
comment = re.compile(r"//! \\ogs_file_(param|attr)\{([A-Za-z_0-9]+)\}( \\todo .*)?$")
comment_special = re.compile(r"//! \\ogs_file(_param|_attr)?_special(\{[A-Za-z_0-9]+\})?( \\todo .*)?$")

# capture #5 is the parameter name
getter = re.compile(r'(get|check|ignore|peek)Config(Parameter|Attribute|Subtree)(List|Optional|All)?'
                   +r'\s*(<.*>)?'
                   +r'\s*\(\s*"([a-zA-Z_0-9:]+)"\s*[,)]')

getter_special = re.compile(r'(get|check|ignore|peek)Config(Parameter|Attribute|Subtree)(List|Optional|All)?'
                           +r'\s*(<.*>)?\(')


# merge lines belonging together from grep -A 2 output.
def merge_lines(it):
    buf = ""
    buf_fn = ""
    buf_lno = 0

    for l in it:
        l = l.strip()
        if (not l) or l == "--":
            # separator line
            if buf_fn:
                yield buf_fn, buf_lno, buf
            buf = ""
            buf_fn = ""
            buf_lno = 0
        else:
            m = re.match("(.*)([:-])([0-9]+)([:-])(.*)$", l)
            assert m
            assert m.group(2) == m.group(4)
            fn = m.group(1)
            lno = int(m.group(3))
            line = m.group(5)
            msg =  fn + m.group(2) + str(lno) + m.group(4) + line

            # remove non-doxygen comments
            line = re.sub('/\*[^!*].*\*/|/\*\*/', '', line)
            line = re.sub("//[^!*].*|//$", "", line, 1)

            if buf_fn:
                if m.group(2) == ":":
                    # new location started, yield the old one
                    yield buf_fn, buf_lno, buf
                    buf = line
                    buf_fn = fn
                    buf_lno = lno
                else:
                    # continuation line
                    assert buf_fn == fn
                    buf += " " + line
            elif m.group(2) == ":":
                buf = line
                buf_fn = fn
                buf_lno = lno
            else:
                # continuation line and empty buffer
                pass

            if buf_fn and (comment.search(line) or comment_special.search(line)):
                # make sure nothing can be appended to doxygen comment lines
                yield buf_fn, buf_lno, buf
                buf = ""
                buf_lno = 0
                buf_fn = ""

    if buf_fn:
        yield buf_fn, buf_lno, buf

    return
    yield


state = "getter"
path = ""
lineno = 0
line = ""
tag_path_comment = ""
param_or_attr_comment = ""

for inline in merge_lines(sys.stdin):
    oldpath = path; oldlineno = lineno; oldline = line
    path, lineno, line = inline

    if path != oldpath: debug(path)

    m = comment.search(line)
    if m:
        if state != "getter":
            write_out("UNNEEDED", oldpath, oldlineno, oldline)
        state = "comment"

        param_or_attr_comment = m.group(1)
        tag_path_comment = m.group(2).replace("__", ".")
        debug(" {0:>5}  //! {1}".format(lineno, tag_path_comment))
        tag_name_comment = tag_path_comment.split(".")[-1]

        continue

    m = comment_special.search(line)
    if m:
        if state != "getter":
            write_out("UNNEEDED", oldpath, oldlineno, oldline)
        state = "comment"
        param_or_attr_comment = "special"

        if m.group(1): # param|attr matched
            # second group must not be empty!
            tag_path_comment = m.group(2).strip("{}").replace("__", ".")
            param = tag_path_comment.split(".")[-1]
            paramtype = ""
            method = ""
            write_out("OK", path, lineno, tag_path_comment, param, paramtype, method)
            state = "getter" # reset state s.t. next time a comment is accepted

        continue

    m = getter.search(line)
    if m:
        param = m.group(5)
        paramtype = m.group(4)[1:-1] if m.group(4) else ""
        method = m.group(1) + "Config" + m.group(2) + (m.group(3) or "")

        if state != "comment" or oldpath != path:
            write_out("NODOC", path, lineno, "NONE", param, paramtype, method)
        else:
            debug(" {0:>5}  {1} {2} ".format(lineno, param, paramtype))

            if param != tag_name_comment:
                debug("error: parameter name from comment and code do not match: "
                        + tag_name_comment + " vs. " + param)
                write_out("NODOC", path, lineno, tag_path_comment, param, paramtype, method)
            elif lineno != oldlineno+1:
                debug("error: the associated comment is not on the line preceding this one."
                        + " line numbers {0} vs. {1}".format(oldlineno, lineno))
                write_out("NODOC", path, lineno, tag_path_comment, param, paramtype, method)
            elif param_or_attr_comment == "param" and m.group(2) != "Parameter" and m.group(2) != "Subtree":
                debug("error: comment says param but code says different.")
                write_out("NODOC", path, lineno, tag_path_comment, param, paramtype, method)
            elif param_or_attr_comment == "attr" and m.group(2) != "Attribute":
                debug("error: comment says attr but code says different.")
                write_out("NODOC", path, lineno, tag_path_comment, param, paramtype, method)
            elif param_or_attr_comment == "special":
                debug("error: comment comments a special line.")
                write_out("NODOC", path, lineno, "UNKNOWN", "UNKNOWN", paramtype, method)
            else:
                write_out("OK", path, lineno, tag_path_comment, param, paramtype, method)

        state = "getter"
        continue

    m = getter_special.search(line)
    if m:
        paramtype = m.group(4)[1:-1] if m.group(4) else ""
        method = m.group(1) + "Config" + m.group(2) + (m.group(3) or "")

        if state != "comment" or oldpath != path:
            write_out("NODOC", path, lineno, "NONE", "UNKNOWN", paramtype, method)
        else:
            if lineno != oldlineno+1:
                debug("error: the associated comment is not on the line preceding this one."
                        + " line numbers {0} vs. {1}".format(oldlineno, lineno))
                write_out("NODOC", path, lineno, "UNKNOWN", "UNKNOWN", paramtype, method)
            elif param_or_attr_comment != "special":
                debug("error: comment does not comment a special line.")
                write_out("NODOC", path, lineno, "UNKNOWN", "UNKNOWN", paramtype, method)
            else:
                write_out("SPECIAL", path, lineno, paramtype, method)

        state = "getter"
        continue

    write_out("WRONGIN", path, lineno, line.strip())
    state = "getter" # reset state in order to avoid warnings
