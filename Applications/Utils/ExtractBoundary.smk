# Usage, e.g.:
#   snakemake -s ExtractBoundary.smk -j 1 --configfile $HOME/code/ogs6/build/buildinfo.yaml
#
# buildinfo.yaml contains variables such as Data_BINARY_DIR

output_path = "FileIO"

import os
os.environ["PATH"] += os.pathsep + os.pathsep.join([config['BIN_DIR']])
workdir: f"{config['Data_BINARY_DIR']}/{output_path}"

# "entry point", otherwise one would had to specify output files as snakemake
# arguments
elem_types = ['tri', 'quad']
rule all:
    input:
        expand("square_10_1x1_{type}_boundary_diff.out", type=elem_types)

include: f"{config['SOURCE_DIR']}/scripts/snakemake/modules/meshes.smk"

rule extract_boundary:
    input:
        rules.generate_square_mesh.output
    output:
        "{mesh_name_prefix}_{size}_{lx}x{ly}_{type}_boundary.vtu"
    shell:
        "ExtractBoundary -i {input} -o {output}"

rule vtkdiff:
    input:
        a = rules.extract_boundary.output,
        b = f"{config['Data_SOURCE_DIR']}/{output_path}/{rules.extract_boundary.output}"
    output:
        "{mesh_name_prefix}_{size}_{lx}x{ly}_{type}_boundary_diff.out"
    params:
        check_mesh = True,
        fields = [
            # second field name can be omitted if identical
            ["bulk_node_ids", 0, 0],
            ["bulk_element_ids", 0, 0],
            ["bulk_face_ids", 0, 0]
        ]
    wrapper:
        f"file://{config['SOURCE_DIR']}/scripts/snakemake/vtkdiff"
