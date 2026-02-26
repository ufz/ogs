(use-modules (guix transformations)
             (guix utils)
             (gnu packages base))

(load "manifest-common.scm")

(define current-dir
  (getcwd))

(define transform1
  (options->transformation (append `((with-source unquote
                                                  (string-append
                                                   "ogs-petsc-mkl="
                                                   current-dir)))
                                   (manifest-eigen-transform-options)
                                   '((with-input . "openmpi=openmpi@4.1.6")))))

(packages->manifest (append (manifest-runtime-packages (transform1 (specification->package
                                                                    "ogs-petsc-mkl")))
                            (list (specification->package "openmpi@4.1.6"))))
