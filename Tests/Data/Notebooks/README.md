# OGS Jupyter notebook tips

Images used in Markdown cells need to be in a subdirectory `images` or `figures`! E.g.:

```markdown
![Schematic view of surfing boundary condition benchmark](figures/surfing_schematic.png)
```

Otherwise the images will not appear on the website.

## Git diff / merge tool

Install [nbdime](https://nbdime.readthedocs.io) and enable [git integration](https://nbdime.readthedocs.io/en/latest/vcs.html):

```bash
nbdime config-git --enable --global
```

## testrunner

```bash
# in virtualenv:
pip install -r requirements.txt
PATH=~/code/ogs/build/release/bin:$PATH python testrunner.py --out _out SimpleMechanics.py
```

In Notebook do checks with `assert False` or raise an exception, e.g. `raise SystemExit()`, on failure.
