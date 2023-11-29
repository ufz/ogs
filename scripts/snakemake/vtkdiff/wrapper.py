"""Snakemake wrapper for vtkdiff."""

__author__ = "Lars Bilke"
__copyright__ = "Copyright 2020, OpenGeoSys Community"
__license__ = "BSD"

from pathlib import Path

from snakemake.shell import shell

# ruff: noqa: F821
output = Path(snakemake.output[0])
if output.exists():
    output.unlink()

if snakemake.params.check_mesh:
    shell("vtkdiff {snakemake.input.a} {snakemake.input.b} -m > {snakemake.output[0]}")

for field in snakemake.params.fields:
    field_a = field[0]
    offset = 0
    if len(field) == 4:
        offset = 1
    field_b = field[0 + offset]
    abs_tol = field[1 + offset]
    rel_tol = field[2 + offset]

    shell(
        """
        vtkdiff {snakemake.input.a} {snakemake.input.b} \
          -a {field_a} -b {field_b} \
          --abs {abs_tol} --rel {rel_tol} 2>&1 | tee -a {snakemake.output[0]}
        """
    )
