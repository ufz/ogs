## Notebook-based tests

Environment defined in this directory's `pyproject.toml` and managed via [`uv`](https://github.com/astral-sh/uv).

3 dependency groups which select the used `ogs` wheel:

- `ogs-local`: Builds ogs locally from source
- `ogs-master` Downloads latest ogs wheel from master branch (requires `--upgrade-package ogs`-parameter)
- `ogs-pypi` Downloads latest ogs release from PyPI

```bash
uv run --group ogs-local pytest

uv run --group ogs-master --upgrade-package ogs pytest
uv run --group ogs-pypi
```

`uv` automatically manages a virtual environment in `./.venv`. It locks the environment according to a `uv.lock`-file. To update the file for all packages run `uv lock --upgrade`. To update only a specific package, e.g. on `ogstools`-release run `uv lock --upgrade-package ogstools`.
