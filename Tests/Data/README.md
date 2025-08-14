## Notebook-based tests

Environment for notebook-based tests is defined in this directory's `pyproject.toml` and managed via [`uv`](https://github.com/astral-sh/uv).

3 dependency groups which select the used `ogs` wheel:

- `ogs-local`: Builds ogs locally from source
- `ogs-master` Downloads latest ogs wheel from master branch (requires `--upgrade-package ogs`-parameter)
- `ogs-pypi` Downloads latest ogs release from PyPI

```bash
uv run --group ogs-local pytest

uv run --group ogs-master --upgrade-package ogs pytest
uv run --group ogs-pypi pytest

uv run --group ogs-local Notebooks/testrunner.py --out /tmp Mechanics/Linear/SimpleMechanics.py
```

`uv` automatically manages a virtual environment in `./.venv`. It locks the environment according to a `uv.lock`-file. To update the file for all packages run `uv lock --upgrade`. To update only a specific package, e.g. on `ogstools`-release run `uv lock --upgrade-package ogstools`.

## Regression and maintenance testing

| | Regression test| Maintenance test |
| ------ | ------ | ------ |
| OGS       | **latest**        | stable (tip of master)       |
| Environment       | stable (pinned via lock-file)       |  **latest** |
| Trigger | Merge-request | scheduled or master-build |
| Role | Developer | Maintainer |

**latest** indicates the "object under test"

Environment may be:

- pip (from PyPI)
- CPM (FindPackage, GIT_TAG / VERSION) - maintenance not implemented
- Compiler version - maintenance not implemented
- ...

When to update a "stable" an environment:

- MUST: Change of underlying stack (e.g. software update of `envinfX` test machines)
- MUST: An MR with a change in the request environment (e.g. a new lib or a new version of an already used lib is needed)
- OPTIONAL: On any other occasion (e.g. OGS release)

Typical workflow to update:

1. Trigger: Maintenance test fails
2. Check latest updates of 3rd party packages and their influence (`uv lock --upgrade`)
3. **MR: Quick fix**: Restrict versions in (open) dependencies, add/open ToDo
4. Maintenance test succeeds (this allows other maintenance issues to be found)
5. Investigation and fix of own software / report to 3rd party
6. **MR: Fix** and reopen dependency
