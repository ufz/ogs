(use-modules (guix transformations)
             (guix utils)
             (gnu packages base))

(load "manifest-common.scm")

(define current-dir
  (getcwd))

(define transform1
  (options->transformation `((with-source unquote
                                          (string-append "ogs-mkl="
                                                         current-dir))
                             (with-commit . "eigen=9000b3767770f6dd0f4cfb12f4e19c71921885a4")
                             (without-tests . "eigen"))))

(packages->manifest
 (manifest-runtime-packages
  (transform1 (specification->package "ogs-mkl"))))
