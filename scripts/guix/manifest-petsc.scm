(use-modules (guix transformations)
             (guix utils)
             (gnu packages base))

(define current-dir (getcwd))

(define transform1
  (options->transformation
    `((with-source . ,(string-append "ogs-petsc=" current-dir))
      (with-commit . "eigen=9000b3767770f6dd0f4cfb12f4e19c71921885a4")
      (without-tests . "eigen")
      (with-configure-flag . "vtk=-DVTK_MODULE_USE_EXTERNAL_VTK_eigen=OFF")
      (with-source . "vtk=https://www.vtk.org/files/release/9.3/VTK-9.3.1.tar.gz")
      (with-input ."openmpi=openmpi@4.1.6"))))

(packages->manifest
 (list
  (transform1 (specification->package "ogs-petsc"))
  (specification->package "coreutils")
  (specification->package "bash"))) ; required for squashfs container image