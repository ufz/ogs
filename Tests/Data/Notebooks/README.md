# OGS Juypter notebook tips

## Git diff / merge tool

Install [nbdime](https://nbdime.readthedocs.io) and enable [git integration](https://nbdime.readthedocs.io/en/latest/vcs.html):

```bash
nbdime config-git --enable --global
```

## testrunner

```bash
# in virtualenv:
pip install -r requirements.txt
PATH=~/code/ogs/build/release/bin:$PATH python testrunner.py --out _out SimpleMechanics.ipynb
```

In Notebook do checks with `assert False` or `raise SystemExit()` on failure.
