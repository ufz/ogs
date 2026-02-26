(define (manifest-runtime-packages package)
  (list package
        (specification->package "vtkdiff")
        (specification->package "which")
        (specification->package "coreutils")
        (specification->package "bash"))) ; required for squashfs container image

(define (manifest-eigen-transform-options)
  '((with-commit . "eigen=9000b3767770f6dd0f4cfb12f4e19c71921885a4") ; or: (with-version . "eigen=5.0.1")
    (without-tests . "eigen")))
