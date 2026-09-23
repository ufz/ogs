# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

from markdown_conversion import remove_duplicate_title_header


def test_removes_first_h1_matching_frontmatter_title():
    markdown = '+++\ntitle = "Example"\n+++\n\n# Example\n\nContent\n'

    assert remove_duplicate_title_header(markdown, "Example") == (
        '+++\ntitle = "Example"\n+++\n\n\nContent\n'
    )


def test_keeps_first_h1_with_a_different_title():
    markdown = "# Notebook heading\n\n# Example\n"

    assert remove_duplicate_title_header(markdown, "Example") == markdown


def test_does_not_treat_h1_in_a_code_block_as_markdown():
    markdown = "```python\n# Example\n```\n\n# Example\n"

    assert remove_duplicate_title_header(markdown, "Example") == (
        "```python\n# Example\n```\n\n"
    )
