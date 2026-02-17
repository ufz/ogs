(use-modules (guix transformations)
(guix utils)
(gnu packages base))

(define current-dir
(getcwd))

(define transform1
(options->transformation `((with-source unquote
                             (string-append "ogs-petsc-mkl="
                                            current-dir))
                (with-commit . "eigen=9000b3767770f6dd0f4cfb12f4e19c71921885a4")
                (without-tests . "eigen"))))

(packages->manifest (list (transform1 (specification->package "ogs-petsc-mkl"))
             (specification->package "coreutils")
             (specification->package "bash")))