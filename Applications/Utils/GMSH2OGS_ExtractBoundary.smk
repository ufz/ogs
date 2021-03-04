# Usage, e.g.:
#   snakemake -s GMSH2OGS_ExtractBoundary.smk -j 1 --configfile $HOME/code/ogs6/build/buildinfo.yaml
#
# buildinfo.yaml contains variables such as Data_BINARY_DIR

import os
os.environ["PATH"] += os.pathsep + os.pathsep.join([config['BIN_DIR']])
workdir: f"{config['Data_BINARY_DIR']}/MeshLib"

#import ipdb
#ipdb.set_trace()

VTKDIFF = f"file://{config['SOURCE_DIR']}/scripts/snakemake/vtkdiff"
SOURCE_DIR = f"{config['Data_SOURCE_DIR']}/MeshLib"

MESH_NAME = "A2"
BULK_MESH = f"{MESH_NAME}.vtu"
GMSH_MESH = f"{MESH_NAME}-gmsh.msh"
INNER_BOUNDARIES = "[12]"
OUTER_BOUNDARIES = "[034567]"

ALL_REFERENCE_FILES = glob_wildcards(f"{SOURCE_DIR}/{{file, {MESH_NAME}.*}}.vtu").file

rule all:
    input:
        expand("{file}.vtkdiff.out", file=ALL_REFERENCE_FILES)

rule gmsh2ogs:
    input:
        f"{SOURCE_DIR}/{GMSH_MESH}"
    output:
        expand("{file}.vtu", file=ALL_REFERENCE_FILES)
    shell:
        "GMSH2OGS -i {input} -o {BULK_MESH} -e -b -v"

rule vtkdiff_geometry_bulk:
    input:
        a = f"{SOURCE_DIR}/{BULK_MESH}",
        b = BULK_MESH
    output:
        f"{MESH_NAME}.vtkdiff.out"
    params:
        check_mesh = True,
        fields = [["MaterialIDs", 0, 0]]
    wrapper:
        VTKDIFF

rule vtkdiff_geometry_boundary:
    input:
        a = f"{SOURCE_DIR}/{MESH_NAME}_{{x}}.vtu",
        b = f"{MESH_NAME}_{{x}}.vtu"
    output:
        f"{MESH_NAME}_{{x}}.vtkdiff.out-bulk_node_ids"
    params:
        check_mesh = True,
        fields = [["bulk_node_ids", 0, 0]]
    wrapper:
        VTKDIFF

rule vtkdiff_geometry_inner:
    input:
        f"{MESH_NAME}_{{x}}.vtkdiff.out-bulk_node_ids",
        a = f"{SOURCE_DIR}/{MESH_NAME}_{{x}}.vtu",
        b = f"{MESH_NAME}_{{x}}.vtu"
    output:
        f"{MESH_NAME}_{{x, {INNER_BOUNDARIES}}}.vtkdiff.out"
    params:
        check_mesh = True,
        fields = [["number_bulk_elements", 0, 0]]
    wrapper:
        VTKDIFF

rule vtkdiff_geometry_outer:
    input:
        f"{MESH_NAME}_{{x}}.vtkdiff.out-bulk_node_ids",
        a = f"{SOURCE_DIR}/{MESH_NAME}_{{x}}.vtu",
        b = f"{MESH_NAME}_{{x}}.vtu"
    output:
        f"{MESH_NAME}_{{x, {OUTER_BOUNDARIES}}}.vtkdiff.out"
    params:
        check_mesh = True,
        fields = [["bulk_element_ids", 0, 0]]
    wrapper:
        VTKDIFF
