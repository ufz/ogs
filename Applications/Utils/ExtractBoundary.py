# Usage, e.g.:
#   python3 ExtractBoundary.py ./buildinfo.yaml
#
# buildinfo.yaml contains variables such as Data_BINARY_DIR

import os, parsl, sys, yaml
from parsl import python_app, bash_app
from parsl.data_provider.files import File

output_path = "FileIO"
elem_types = ["tri", "quad"]

parsl.load()

config = dict()
with open(sys.argv[1]) as f:
    config = yaml.safe_load(f)

os.environ["PATH"] += os.pathsep + os.pathsep.join([config["BIN_DIR"]])
os.chdir(f"{config['Data_BINARY_DIR']}/{output_path}")

print(f"{config['Data_BINARY_DIR']}/{output_path}")


# Apps
@bash_app
def generate_meshes(
    elem_type, outputs=[], stderr=parsl.AUTO_LOGNAME, stdout=parsl.AUTO_LOGNAME
):
    return f"""generateStructuredMesh -e {elem_type} \
            --lx 1 --ly 1 \
            --nx 10 --ny 10 \
            -o input_square_1x1_{elem_type}.vtu"""


@bash_app
def extract_boundary(
    elem_type, inputs=[], outputs=[], stderr=parsl.AUTO_LOGNAME, stdout=parsl.AUTO_LOGNAME
):
    return f"""ExtractBoundary -i {inputs[0].filepath} \
            -o square_1x1_{elem_type}_boundary.vtu"""


# compares the files in inputs[0] and inputs[1]
@bash_app
def vtk_diff(
    fields, inputs=[], outputs=[], stderr=parsl.AUTO_LOGNAME, stdout=parsl.AUTO_LOGNAME
):
    import os

    script = ""
    if os.path.exists(outputs[0]):
        os.remove(outputs[0])
    for field in fields:
        field_a = field[0]
        offset = 0
        if len(field) == 4:
            offset = 1
        field_b = field[0 + offset]
        abs_tol = field[1 + offset]
        rel_tol = field[2 + offset]

        script += f"""vtkdiff {inputs[0]} {inputs[1]} \
              -a {field_a} -b {field_b} \
              --abs {abs_tol} --rel {rel_tol} >> {outputs[0]}
        """
    return script


# Workflow
for elem_type in elem_types:
    gm = generate_meshes(elem_type, outputs=[File(f"input_square_1x1_{elem_type}.vtu")])
    eb = extract_boundary(
        elem_type,
        inputs=[gm.outputs[0]],
        outputs=[File(f"square_1x1_{elem_type}_boundary.vtu")],
    )
    diff = vtk_diff(
        fields=[
            # second field name can be omitted if identical
            ["bulk_node_ids", 0, 0],
            ["bulk_element_ids", 0, 0],
            ["bulk_face_ids", 0, 0],
        ],
        inputs=[
            eb.outputs[0],
            f"{config['Data_SOURCE_DIR']}/{output_path}/{eb.outputs[0].filename}",
        ],
        outputs=[File(f"square_1x1_{elem_type}_boundary_diff.out")],
    )
    print(diff.result())
