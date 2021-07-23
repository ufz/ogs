+++
date = "2018-02-26T11:00:13+01:00"
title = "Code style and formatting"
author = "Lars Bilke"
weight = 1014

[menu]
  [menu.devguide]
    parent = "development-workflows"
+++

We aim for a consistent and readable coding style. You do not need to worry about styling if you use the right tools we present in the following. Please also enable [pre-commit]({{< ref "setup-fork.md#optional-enable-git-commit-hooks" >}}) in your repository to have these check run automatically for you on every commit!

## C++

Use [`clang-format`](https://clang.llvm.org/docs/ClangFormat.html). It can be added to be automatically run on save in your editor / IDE:

- [Vim](https://github.com/rhysd/vim-clang-format)
- [Visual Studio](https://devblogs.microsoft.com/cppblog/clangformat-support-in-visual-studio-2017-15-7-preview-1/)
- [Visual Studio Code](https://marketplace.visualstudio.com/items?itemName=xaver.clang-format)

The current style is defined in [.clang-format](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/.clang-format).

## Python

Use [`black`](https://black.readthedocs.io/en/stable/). It can be added to be automatically run on save in your editor / IDE:

- [Vim](https://black.readthedocs.io/en/stable/editor_integration.html#vim)
- [PyCharm](https://black.readthedocs.io/en/stable/editor_integration.html#pycharm-intellij-idea)
- [Visual Studio Code](https://code.visualstudio.com/docs/python/editing#_formatting)

`black` is also run by our `pre-commit`-hooks. To run manually:

```bash
pre-commit run black --all-files
```

## CMake

Use [cmake-format](https://cmake-format.readthedocs.io/en/latest/cmake-format.html). It can be added to be automatically run on save in your editor / IDE:

- [Visual Studio Code](https://marketplace.visualstudio.com/items?itemName=cheshirekow.cmake-format)

The current style is defined in [.cmake-format.yaml](https://gitlab.opengeosys.org/ogs/ogs/-/blob/master/.cmake-format.yaml).

On installation make sure to install with the yaml option:

```
pip install cmakelang[YAML]
```

## Build targets

You can run `clang-format` and `cmake-format` on the OGS code base by building the following targets:

| Language | Shows formatting | Applies formatting |
| -------- | ---------------- | ------------------ |
| C++      | `clang-format`   | `fix-clang-format` |
| CMake    | `cmake-format`   | `fix-cmake-format` |
