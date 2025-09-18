(use-modules (guix transformations)
(guix utils)
(gnu packages base))

(define current-dir
(getcwd))

(define transform1
(options->transformation `((with-source unquote
                             (string-append "ogs-serial="
                                            current-dir))
                (with-commit . "eigen=9000b3767770f6dd0f4cfb12f4e19c71921885a4")
                (without-tests . "eigen")
                (with-configure-flag . "vtk=-DVTK_MODULE_USE_EXTERNAL_VTK_eigen=OFF")
                (with-commit . "vtk=v9.3.1")
                (with-git-url . "vtk=https://gitlab.kitware.com/vtk/vtk.git"))))

(packages->manifest (list (transform1 (specification->package "ogs-petsc-mkl"))
             (specification->package "coreutils")
             (specification->package "bash")))