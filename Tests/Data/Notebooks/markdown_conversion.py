# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

import re

_CODE_FENCE = re.compile(r"^ {0,3}(`{3,}|~{3,})")
_H1_HEADER = re.compile(r"^ {0,3}#(?!#)[ \t]+(.*?)[ \t]*\r?\n?$")


def remove_duplicate_title_header(markdown, title):
    """Remove the first Markdown H1 if it duplicates the frontmatter title."""
    if not title:
        return markdown

    lines = markdown.splitlines(keepends=True)
    code_fence = None
    for line_number, line in enumerate(lines):
        fence_match = _CODE_FENCE.match(line)
        if code_fence is not None:
            if (
                fence_match is not None
                and fence_match.group(1)[0] == code_fence[0]
                and len(fence_match.group(1)) >= code_fence[1]
            ):
                code_fence = None
            continue

        if fence_match is not None:
            marker = fence_match.group(1)
            code_fence = (marker[0], len(marker))
            continue

        heading_match = _H1_HEADER.match(line)
        if heading_match is None:
            continue

        if heading_match.group(1).strip() == title:
            del lines[line_number]
        break

    return "".join(lines)
