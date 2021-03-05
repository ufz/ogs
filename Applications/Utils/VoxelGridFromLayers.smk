# Usage, e.g.:
#   snakemake -s VoxelGridFromLayers.smk -j 1 --configfile $HOME/code/ogs6/build/buildinfo.yaml
#
# buildinfo.yaml contains variables such as Data_BINARY_DIR

import os
os.environ["PATH"] += os.pathsep + os.pathsep.join([config['BIN_DIR']])
workdir: f"{config['Data_SOURCE_DIR']}/MeshLib"
out_dir = f"{config['Data_BINARY_DIR']}/MeshLib"

rule all:
    input:
        f"{out_dir}/AREHS_test_diff_geometry.out",
        f"{out_dir}/AREHS_test_fault_diff.out",
        f"{out_dir}/AREHS_test_iso_diff_geometry.out"

rule layers_to_grid:
    input:
        "AREHS_test_layers.txt"
    output:
        f"{out_dir}/AREHS_test.vtu"
    shell:
        "Layers2Grid -i {input} -o {output} -x 500 -y 300 -z 100"

rule vtkdiff_grid_geometry:
    input:
        out = rules.layers_to_grid.output,
        ref = "AREHS_test.vtu"
    output:
        f"{out_dir}/AREHS_test_diff_geometry.out"
    shell:
        "vtkdiff -m {input.out} {input.ref} 2>&1 | tee {output}"

rule add_fault_to_grid:
    input:
        grid = rules.layers_to_grid.output,
        fault = "AREHS_fault.vtu"
    output:
        f"{out_dir}/AREHS_test_fault.vtu"
    shell:
        "AddFaultToVoxelGrid -i {input.grid} -f {input.fault} -o {output}"

rule vtkdiff_fault:
    input:
        a = rules.add_fault_to_grid.output,
        b = "AREHS_test_fault.vtu"
    output:
        f"{out_dir}/AREHS_test_fault_diff.out"
    params:
        check_mesh = True,
        fields = [["MaterialIDs", 0, 0]]
    wrapper:
        f"file://{config['SOURCE_DIR']}/scripts/snakemake/vtkdiff"

rule layers_to_grid_iso:
    input:
        "AREHS_test_layers.txt"
    output:
        f"{out_dir}/AREHS_test_iso.vtu"
    shell:
        "Layers2Grid -i {input} -o {output} -x 500"

rule vtkdiff_grid_iso_geometry:
    input:
        out = rules.layers_to_grid_iso.output,
        ref = "AREHS_test_iso.vtu"
    output:
        f"{out_dir}/AREHS_test_iso_diff_geometry.out"
    shell:
        "vtkdiff -m {input.out} {input.ref} 2>&1 | tee {output}"
