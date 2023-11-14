(define-module (ogs-package)
  #:use-module (guix)
  #:use-module (guix packages)
  #:use-module (guix build-system cmake)
  #:use-module (guix git-download)
  #:use-module ((guix licenses)
                #:prefix license:)
  #:use-module (gnu packages algebra)
  #:use-module (gnu packages boost)
  #:use-module (gnu packages check)
  #:use-module (gnu packages cmake)
  #:use-module (gnu packages compression)
  #:use-module (gnu packages cpp)
  #:use-module (gnu packages image-processing)
  #:use-module (gnu packages logging)
  #:use-module (gnu packages maths)
  #:use-module (gnu packages mpi)
  #:use-module (gnu packages ninja)
  #:use-module (gnu packages pretty-print)
  #:use-module (gnu packages python-xyz)
  #:use-module (gnu packages python)
  #:use-module (gnu packages version-control)
  #:use-module (gnu packages xml)
  #:use-module (ogs-dependencies))

(define-public ogs
  (package
    (name "ogs")
    (version "6.4.4")
    (source (origin
              (method git-fetch)
              (uri (git-reference
                    (url "https://gitlab.opengeosys.org/ogs/ogs.git")
                    (commit "d4ca7e627f2fc012bfe434649b797e78e5c2a8f1")))
              (file-name (git-file-name name version))
              (sha256
               (base32
                "1f6mcbjx76irf1g0xkh6hgpv4qn2swbiyvlazvlrhjfyxb9bckq9"))))
    (build-system cmake-build-system)
    (arguments
     `(#:build-type "Release"
       #:configure-flags (list) ;empty list to be appended in inherited packages
       #:cmake ,cmake)) ;for newer CMake version
    (inputs (list boost
                  eigen
                  exprtk
                  hdf5
                  iphreeqc
                  json-modern-cxx
                  libxml2
                  metis
                  pybind11-2.10.4
                  python
                  range-v3
                  spdlog
                  tclap
                  tetgen
                  zlib
                  vtk
                  xmlpatch))
    (native-inputs (list git ninja))
    (synopsis "OpenGeoSys")
    (description
     "Simulation of thermo-hydro-mechanical-chemical (THMC) processes in porous and fractured media")
    (home-page "https://www.opengeosys.org")
    (properties '((tunable? . #t)))
    (license license:bsd-3)))

(define-public ogs-ssd
  (package
    (inherit ogs)
    (name "ogs-ssd")
    (arguments
     (substitute-keyword-arguments (package-arguments ogs)
       ((#:configure-flags flags)
        `(cons* "-DOGS_BUILD_PROCESSES=SteadyStateDiffusion"
                ,flags))))
    (synopsis "OGS with SteadyStateDiffusion only (for faster build testing)")))

(define-public ogs-petsc
  (package
    (inherit ogs)
    (name "ogs-petsc")
    (inputs (modify-inputs (package-inputs ogs)
              (prepend openmpi petsc-openmpi)
              (replace "vtk" vtk-openmpi)
              (replace "hdf5" hdf5-parallel-openmpi)))
    (arguments
     (substitute-keyword-arguments (package-arguments ogs)
       ((#:configure-flags flags)
        `(cons* "-DOGS_USE_PETSC=ON" "-DCMAKE_C_COMPILER=mpicc"
                "-DCMAKE_CXX_COMPILER=mpic++"
                ,flags))))
    (synopsis "OGS with PETSc")))

(define-public ogs-petsc-ssd
  (package
    (inherit ogs-petsc)
    (name "ogs-petsc-ssd")
    (arguments
     (substitute-keyword-arguments (package-arguments ogs-petsc)
       ((#:configure-flags flags)
        `(cons* "-DOGS_BUILD_PROCESSES=SteadyStateDiffusion"
                ,flags))))
    (synopsis
     "OGS with PETSc and SteadyStateDiffusion only (for faster build testing)")))

;; return package
ogs
