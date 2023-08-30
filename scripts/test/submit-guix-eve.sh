#!/bin/bash
module load GCC/12.2.0 OpenMPI/4.1.4
echo ${prj_file}
srun apptainer exec /data/ogs/apptainer/guix/ogs-petsc-ssd_head.squashfs ogs ${prj_file}
