(use-modules (guix gexp)
             (guix packages)
             (guix utils)
             (gnu packages maths))

(define gcc-toolchain-15
  (specification->package "gcc-toolchain@15.2.0"))

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

(define (package-with-gcc-toolchain-15 pkg)
  (package
    (inherit pkg)
    (native-inputs (modify-inputs (package-native-inputs pkg)
                    (append gcc-toolchain-15)))))

(define (manifest-runtime-packages package)
  (list (package-with-metis-int64 (package-with-gcc-toolchain-15 package))
        (specification->package "python-wrapper")
        (specification->package "vtkdiff")
        (specification->package "which")
        (specification->package "coreutils")
        (specification->package "bash")))

(define (manifest-eigen-transform-options)
  ; OpenMP fix by Jörg
  ; https://gitlab.com/libeigen/eigen/-/merge_requests/1786
  '((with-commit . "eigen=24e0c2a125d2b37b35719124d1f758777c150ca8")
    (without-tests . "eigen")))
