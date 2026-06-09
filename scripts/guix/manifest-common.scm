(use-modules (guix gexp)
             (guix packages)
             (guix utils)
             (gnu packages maths))

(define metis-int64
  (package
    (inherit metis)
    (arguments
     (substitute-keyword-arguments (package-arguments metis)
       ((#:phases phases
         #~%standard-phases)
        #~(modify-phases #$phases
            (add-after 'unpack 'set-idx-width
              (lambda _
                (substitute* "include/metis.h"
                  (("#define IDXTYPEWIDTH 32")
                   "#define IDXTYPEWIDTH 64"))))))))))

(define (package-with-metis-int64 pkg)
  (package
    (inherit pkg)
    (inputs (modify-inputs (package-inputs pkg)
              (replace "metis" metis-int64)))
    (propagated-inputs (modify-inputs (package-propagated-inputs pkg)
                         (replace "metis" metis-int64)))))

(define (manifest-runtime-packages package)
  (list (package-with-metis-int64 package)
        (specification->package "python-wrapper")
        (specification->package "vtkdiff")
        (specification->package "which")
        (specification->package "coreutils")
        (specification->package "bash")))

(define (manifest-eigen-transform-options)
  '((with-commit . "eigen=9000b3767770f6dd0f4cfb12f4e19c71921885a4")
    (without-tests . "eigen")))
