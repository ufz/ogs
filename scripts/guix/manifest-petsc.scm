(use-modules (guix transformations)
             (guix utils)
             (gnu packages base))

(define current-dir (getcwd))

(define transform1
  (options->transformation
    `((with-source . ,(string-append "ogs-petsc=" current-dir))
      (with-commit . "eigen=9000b3767770f6dd0f4cfb12f4e19c71921885a4")
      (without-tests . "eigen")
      (with-input ."openmpi=openmpi@4.1.6"))))

(packages->manifest
 (list
  (transform1 (specification->package "ogs-petsc"))
  (specification->package "coreutils")
  (specification->package "bash"))) ; required for squashfs container image