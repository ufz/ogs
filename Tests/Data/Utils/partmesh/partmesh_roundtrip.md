+++

+++

## `partmesh` round-trip

```python
import os
from pathlib import Path

out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
out_dir.mkdir(parents=True, exist_ok=True)

data_dir = os.environ.get('OGS_DATA_DIR')

input_mesh_basename = "cube_1x1x1_hex_1e3"
input_mesh = f"{data_dir}/EllipticPETSc/{input_mesh_basename}.vtu"
num_partitions = 4
```

```python
! partmesh --ogs2metis -i {input_mesh} -o {out_dir}
```

```python
! partmesh --exe_metis -n {num_partitions} -i {input_mesh} -o {out_dir}
```

```python
! cd {out_dir} && mpirun --bind-to none -np {num_partitions} binaryToPVTU -i {input_mesh_basename} -o {input_mesh_basename}
```

Please note that `binaryToPVTU` has to be run with MPI and is therefore available on OGS PETSc configurations only.
See also https://www.opengeosys.org/docs/tools/meshing/reordermesh/#why-is-the-tool-necessary for an explanation of the following.

```python
! cd {out_dir} && pvtu2vtu -i {input_mesh_basename}.pvtu -o {input_mesh_basename}.vtu
! cd {out_dir} && identifySubdomains -m {input_mesh} -o identify_ -- {input_mesh_basename}.vtu
! ReorderMesh -i {out_dir}/identify_{input_mesh_basename}.vtu -o {out_dir}/{input_mesh_basename}_reordered.vtu
! vtkdiff --mesh_check -- {input_mesh} {out_dir}/{input_mesh_basename}_reordered.vtu
assert _exit_code == 0 # last exit code (from vtkdiff call)
```
