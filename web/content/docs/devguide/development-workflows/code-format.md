+++
date = "2018-02-26T11:00:13+01:00"
title = "Code style and formatting"
author = "Lars Bilke"
weight = 1012

[menu]
  [menu.devguide]
    parent = "development-workflows"
+++

We aim for a consistent and readable coding style. You do not need to worry about styling if you use the right tools we present in the following.

## C++

Use [`clang-format`](https://clang.llvm.org/docs/ClangFormat.html). Can be added to be automatically run on save in your editor / IDE:

- [Vim](https://github.com/rhysd/vim-clang-format)
- [Visual Studio](https://devblogs.microsoft.com/cppblog/clangformat-support-in-visual-studio-2017-15-7-preview-1/)
- [Visual Studio Code](https://marketplace.visualstudio.com/items?itemName=xaver.clang-format)

Current style is defined in [.clang-format](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/.clang-format).

A pre-commit hook for git checking the code formatting can be found [here](https://gitlab.opengeosys.org/ogs/ogs-utils/-/tree/master/dev/code-formatting/clang-format-pre-commit-hook).

## Python

Use [`black`](https://black.readthedocs.io/en/stable/). Can be added to be automatically run on save in your editor / IDE:

- [Vim](https://black.readthedocs.io/en/stable/editor_integration.html#vim)
- [PyCharm](https://black.readthedocs.io/en/stable/editor_integration.html#pycharm-intellij-idea)
- [Visual Studio Code](https://code.visualstudio.com/docs/python/editing#_formatting)

`black` is also run by our `pre-commit`-hooks. To run manually:

```bash
pre-commit run black --all-files
```
