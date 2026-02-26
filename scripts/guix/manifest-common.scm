(define (manifest-runtime-packages package)
  (list package
        (specification->package "vtkdiff")
        (specification->package "which")
        (specification->package "coreutils")
        (specification->package "bash"))) ; required for squashfs container image
