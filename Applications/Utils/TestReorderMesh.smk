# buildinfo.yaml contains variables such as Data_BINARY_DIR

output_path = "FileIO"

import os
os.environ["PATH"] = os.pathsep.join([config['BIN_DIR']]) + os.pathsep + os.environ["PATH"]
workdir: f"{config['Data_BINARY_DIR']}/{output_path}"

input_data_dir = f"{config['SOURCE_DIR']}/Tests/Data/EllipticPETSc/ReorderingInputData/"

refinements = [2,4]
element_types = ['hex', 'prism', 'tet', 'pyramid']
number_of_partitions = [2]
prj_base_name = 'steady_state'

rule all:
    input:
        expand("cube_1x1x1_{type}_{refine}x{refine}x{refine}_{number_of_partitions}_vtk_diff.out", type=element_types, refine=refinements, number_of_partitions=number_of_partitions)

rule generate_hex_mesh:
    output:
        "cube_1x1x1_{type}_{lx,\d+}x{ly,\d+}x{lz,\d+}.vtu"
    shell:
        """
        generateStructuredMesh -e {wildcards.type} \
            --lx 1 --ly 1 --lz 1 \
            --nx {wildcards.lx} --ny {wildcards.ly}  --nz {wildcards.lz} \
            -o {output}
        """

rule partitioning:
    input:
        a = rules.generate_hex_mesh.output
    output:
        "cube_1x1x1_{type}_{lx}x{ly}x{lz}_{number_of_partitions}/cube_1x1x1_{type}_{lx}x{ly}x{lz}.mesh"
    shell:
        """
        partmesh -s -i {input.a} -o cube_1x1x1_{wildcards.type}_{wildcards.lx}x{wildcards.ly}x{wildcards.lz}_{wildcards.number_of_partitions}

        partmesh -m -n {wildcards.number_of_partitions} -i {input.a} -o cube_1x1x1_{wildcards.type}_{wildcards.lx}x{wildcards.ly}x{wildcards.lz}_{wildcards.number_of_partitions}
        cd cube_1x1x1_{wildcards.type}_{wildcards.lx}x{wildcards.ly}x{wildcards.lz}_{wildcards.number_of_partitions}

        mv cube_1x1x1_{wildcards.type}_{wildcards.lx}x{wildcards.ly}x{wildcards.lz}_partitioned_cell_properties_cfg{wildcards.number_of_partitions}.bin bulk_mesh_partitioned_cell_properties_cfg{wildcards.number_of_partitions}.bin
        mv cube_1x1x1_{wildcards.type}_{wildcards.lx}x{wildcards.ly}x{wildcards.lz}_partitioned_cell_properties_val{wildcards.number_of_partitions}.bin bulk_mesh_partitioned_cell_properties_val{wildcards.number_of_partitions}.bin
        mv cube_1x1x1_{wildcards.type}_{wildcards.lx}x{wildcards.ly}x{wildcards.lz}_partitioned_msh_cfg{wildcards.number_of_partitions}.bin bulk_mesh_partitioned_msh_cfg{wildcards.number_of_partitions}.bin
        mv cube_1x1x1_{wildcards.type}_{wildcards.lx}x{wildcards.ly}x{wildcards.lz}_partitioned_msh_ele{wildcards.number_of_partitions}.bin bulk_mesh_partitioned_msh_ele{wildcards.number_of_partitions}.bin
        mv cube_1x1x1_{wildcards.type}_{wildcards.lx}x{wildcards.ly}x{wildcards.lz}_partitioned_msh_ele_g{wildcards.number_of_partitions}.bin bulk_mesh_partitioned_msh_ele_g{wildcards.number_of_partitions}.bin
        mv cube_1x1x1_{wildcards.type}_{wildcards.lx}x{wildcards.ly}x{wildcards.lz}_partitioned_msh_nod{wildcards.number_of_partitions}.bin bulk_mesh_partitioned_msh_nod{wildcards.number_of_partitions}.bin
        """

rule copy_project_and_gml_files:
    input:
        project_file_in = input_data_dir + prj_base_name + ".prj",
        cube_output_xml_file_in = input_data_dir + "cube_output.xml",
        steady_state_diffusion_xml_file_in = input_data_dir + "SteadyStateDiffusion.xml",
        gml_file_in = input_data_dir + "cube_1x1x1.gml"
    output:
        project_file_out = "cube_1x1x1_{type}_{lx}x{ly}x{lz}_{number_of_partitions}/" + prj_base_name + ".prj",
        cube_output_xml_file_out = "cube_1x1x1_{type}_{lx}x{ly}x{lz}_{number_of_partitions}/cube_output.xml",
        steady_state_diffusion_xml_file_out = "cube_1x1x1_{type}_{lx}x{ly}x{lz}_{number_of_partitions}/SteadyStateDiffusion.xml",
        gml_file_out = "cube_1x1x1_{type}_{lx}x{ly}x{lz}_{number_of_partitions}/cube_1x1x1.gml"
    shell:
        """
        cp {input.project_file_in} {output.project_file_out}
        cp {input.cube_output_xml_file_in} {output.cube_output_xml_file_out}
        cp {input.steady_state_diffusion_xml_file_in} {output.steady_state_diffusion_xml_file_out}
        cp {input.gml_file_in} {output.gml_file_out}
        """

rule execute_ogs:
    input:
        a = rules.copy_project_and_gml_files.output,
        b = rules.partitioning.output
    output:
        "cube_1x1x1_{type}_{lx}x{ly}x{lz}_{number_of_partitions}/results/bulk_mesh_ts_1_t_1_000000.pvtu"
    shell:
        """
        (
        cd cube_1x1x1_{wildcards.type}_{wildcards.lx}x{wildcards.ly}x{wildcards.lz}_{wildcards.number_of_partitions}
        mpirun -np {wildcards.number_of_partitions} ogs {prj_base_name}.prj -o results/
        )
        """

rule remove_ghost_data:
    input:
        rules.execute_ogs.output
    output:
        "cube_1x1x1_{type}_{lx}x{ly}x{lz}_{number_of_partitions}/results/bulk_mesh_ts_1_t_1_000000.vtu"
    shell:
        """
        pvtu2vtu -i {input} -o {output}
        """

rule identifySubdomains:
    input:
        orig = rules.generate_hex_mesh.output,
        modified = rules.remove_ghost_data.output
    output:
        "cube_1x1x1_{type}_{lx}x{ly}x{lz}_{number_of_partitions}/results/with_ids_bulk_mesh_ts_1_t_1_000000.vtu"
    shell:
        """
        cd cube_1x1x1_{wildcards.type}_{wildcards.lx}x{wildcards.ly}x{wildcards.lz}_{wildcards.number_of_partitions}/results
        identifySubdomains -m ../../{input.orig} -f -o with_ids_ -- bulk_mesh_ts_1_t_1_000000.vtu
        """

rule reorder_nodes_elements:
    input:
        rules.identifySubdomains.output
    output:
        "cube_1x1x1_{type}_{lx}x{ly}x{lz}_{number_of_partitions}/results/bulk_mesh_reordered.vtu"
    shell:
        """
        ReorderMesh -i {input} -o {output}
        """

rule vtk_diff:
    input:
        orig = rules.generate_hex_mesh.output,
        reordered = rules.reorder_nodes_elements.output
    output:
        "cube_1x1x1_{type}_{lx}x{ly}x{lz}_{number_of_partitions}_vtk_diff.out"
    shell:
        """
        vtkdiff -m {input.orig} {input.reordered} > {output}
        """
