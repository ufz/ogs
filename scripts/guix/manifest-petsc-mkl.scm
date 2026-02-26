(use-modules (guix transformations)
(guix utils)
(gnu packages base))

(load "manifest-common.scm")

(define current-dir
(getcwd))

(define transform1
  (options->transformation
    (append
     `((with-source . ,(string-append "ogs-petsc-mkl=" current-dir)))
     (manifest-eigen-transform-options))))

(packages->manifest
 (manifest-runtime-packages
  (transform1 (specification->package "ogs-petsc-mkl"))))
