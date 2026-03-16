;; Will be added upstream after 6.5.8 release.
(define-module (ogs-xdmf)
  #:use-module (guix)
  #:use-module (guix packages)
  #:use-module (guix download)
  #:use-module (guix git-download)
  #:use-module (guix gexp)
  #:use-module (guix build-system cmake)
  #:use-module (ogs-mkl)
  #:use-module (gnu packages geo)
  #:use-module (gnu packages boost)
  #:use-module (gnu packages xml)
  #:use-module (gnu packages maths)
  #:use-module ((guix licenses) #:prefix license:))

(define-public xdmf
  (package
    (name "xdmf")
    (version (git-version "3.0.0" "1" "04a84bab0eb1568e0f1a27c8fb60c6931efda003"))
    (source
      (origin
        (method git-fetch)
        (uri (git-reference
               (url "https://gitlab.kitware.com/xdmf/xdmf")
               (commit "04a84bab0eb1568e0f1a27c8fb60c6931efda003")))
        (file-name (git-file-name name version))
        (patches
          (list
            (origin
              (method url-fetch)
              ;; Patch for newer HDF5 versions
              (uri "https://gitlab.opengeosys.org/ogs/xdmflib/-/commit/92a851f1acb87ad5367eb62f9b97785bedb700bb.patch")
              (sha256
                (base32 "0gwh42saxww00j3432b7pbimb3dzzzayha2pnvvjn2insbbx835c")))))
        (sha256
          (base32 "06k4vibkvgxlzkn06x470aq5q18p7yhql8awrpdz3czys3z8i288"))))
    (build-system cmake-build-system)
    (inputs
      (list boost hdf5 libxml2))
    (arguments
      (list
        #:tests? #f)) ;; has no tests
    (home-page "https://www.xdmf.org/index.html")
    (synopsis "XDMF library")
    (description
      "This package provides the eXtensible Data Model and Format (XDMF) C++ library.")
    (license license:bsd-4)))

(define-public ogs-xdmf
  (package
    (inherit ogs-serial)
    (name "ogs-xdmf")
    (inputs
      (modify-inputs (package-inputs ogs-serial)
        (prepend xdmf)))
    (synopsis "OpenGeoSys with xdmfdiff")))

(define-public ogs-petsc-xdmf
  (package
    (inherit ogs-petsc)
    (name "ogs-petsc-xdmf")
    (inputs
      (modify-inputs (package-inputs ogs-petsc)
        (prepend xdmf)))
    (synopsis "OpenGeoSys with PETSc and xdmfdiff")))

(define-public ogs-mkl-xdmf
  (package
    (inherit ogs-serial)
    (name "ogs-mkl-xdmf")
    (inputs
      (modify-inputs (package-inputs ogs-mkl)
        (prepend xdmf)))
    (synopsis "OpenGeoSys with MKL and xdmfdiff")))

(define-public ogs-petsc-mkl-xdmf
  (package
    (inherit ogs-petsc)
    (name "ogs-petsc-mkl-xdmf")
    (inputs
      (modify-inputs (package-inputs ogs-petsc-mkl)
        (prepend xdmf)))
    (synopsis "OpenGeoSys with PETSc, MKL, and xdmfdiff")))
