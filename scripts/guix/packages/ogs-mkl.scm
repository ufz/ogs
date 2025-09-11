(define-module (ogs-mkl)
  #:use-module (guix)
  #:use-module (guix packages)
  #:use-module (gnu packages geo)
  #:use-module (guix-science-nonfree packages mkl)
  #:use-module (guix-science-nonfree packages maths))

(define-public ogs-mkl
  (package
    (inherit ogs-serial)
    (name "ogs-mkl")
    (inputs (modify-inputs (package-inputs ogs-serial)
              (prepend intel-oneapi-mkl)))
    (arguments
     (substitute-keyword-arguments (package-arguments ogs-serial)
       ((#:configure-flags flags)
        #~(cons* "-DOGS_USE_MKL=ON"
                 (string-append "-DMKLROOT="
                                (assoc-ref %build-inputs "intel-oneapi-mkl"))
                 #$flags))))
    (synopsis "OGS with MKL (sequential!)")))

(define-public ogs-petsc-mkl
  (package
    (inherit ogs-petsc)
    (name "ogs-petsc-mkl")
    (inputs (modify-inputs (package-inputs ogs-petsc)
              (prepend intel-oneapi-mkl)
              (replace "petsc-openmpi" petsc-openmpi-mkl)))
    (arguments
     (substitute-keyword-arguments (package-arguments ogs-petsc)
       ((#:configure-flags flags)
        #~(cons* "-DOGS_USE_MKL=ON"
                 (string-append "-DMKLROOT="
                                (assoc-ref %build-inputs "intel-oneapi-mkl"))
                 #$flags))))
    (synopsis "OGS with PETSc and MKL (sequential only)")))
