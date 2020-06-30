# Usage, e.g.:
#   # generateStructuredMesh and ExtractBoundary have to be in the path
#   snakemake -s ExtractBoundary.smk -j 1 --configfile $HOME/code/ogs6/build/Tests/snakemake.yaml -d $HOME/code/ogs6/build/Tests/Data/FileIO

import os
os.environ["PATH"] += os.pathsep + os.pathsep.join([config['BIN_DIR']])

# "entry point", otherwise one would had to specify output files as snakemake
# arguments
elem_types = ['tri', 'quad']
rule all:
    input:
        expand("square_1x1_{type}_boundary_diff.out", type=elem_types)

rule generate_meshes:
    output:
        "input_square_1x1_{type}.vtu"
    shell:
        """
        generateStructuredMesh -e {wildcards.type} \
            --lx 1 --ly 1 \
            --nx 10 --ny 10 \
            -o {output}
        """

rule extract_boundary:
    input:
        "input_square_1x1_{type}.vtu"
    output:
        "square_1x1_{type}_boundary.vtu"
    shell:
        "ExtractBoundary -i {input} -o {output}"

rule vtkdiff:
    input:
        "square_1x1_{type}_boundary.vtu"
    output:
        "square_1x1_{type}_boundary_diff.out"
    shell:
        """
        vtkdiff {input} {config[Data_SOURCE_DIR]}/FileIO/square_1x1_{wildcards.type}_boundary.vtu -a bulk_node_ids -b bulk_node_ids --abs 0 --rel 0 > {output}
        vtkdiff {input} {config[Data_SOURCE_DIR]}/FileIO/square_1x1_{wildcards.type}_boundary.vtu -a bulk_element_ids -b bulk_element_ids --abs 0 --rel 0 >> {output}
        vtkdiff {input} {config[Data_SOURCE_DIR]}/FileIO/square_1x1_{wildcards.type}_boundary.vtu -a bulk_face_ids -b bulk_face_ids --abs 0 --rel 0 >> {output}
        """
