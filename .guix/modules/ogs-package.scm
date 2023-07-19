;; This file defines a Guix package.  It can be used to spawn an
;; interactive development environment:
;;
;;   guix shell
;;
;; Or it can be used to build Guile from a checkout in an isolated
;; environment:
;;
;;   guix build -f guix.scm
;;
;; Likewise, you may cross-compile it:
;;
;;   guix build -f guix.scm --target=x86_64-w64-mingw32
;;
;; Or may want to build a variant:
;;   guix build -L $PWD/.guix/modules ogs-petsc-ssd

(define-module (ogs-package)
  #:use-module (guix)
  #:use-module (guix packages)
  #:use-module (guix build-system cmake)
  #:use-module (guix download)
  #:use-module (guix git-download)
  #:use-module ((guix licenses) #:prefix license:)
  #:use-module (gnu packages algebra)
  #:use-module (gnu packages boost)
  #:use-module (gnu packages check)
  #:use-module (gnu packages certs) ; TODO: cpm
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
  )

(define vcs-file?
  ;; Return true if the given file is under version control.
  (or (git-predicate "../..") ; (current-source-directory)
      (const #t)))                                ;not in a Git checkout

(define-public ogs
  (package
    (name "ogs")
    ; (version "6.4.4")
    ; (source (origin
    ;           (method git-fetch)
    ;           (uri (git-reference
    ;                 (url "https://gitlab.opengeosys.org/ogs/ogs.git")
    ;                 (commit "d4ca7e627f2fc012bfe434649b797e78e5c2a8f1")
    ;                 (recursive? #t)))
    ;           (file-name (git-file-name name version))
    ;           (sha256
    ;            ;; Update with `guix hash -rx .`, make sure to have submodules updated!
    ;            (base32
    ;             "1f6mcbjx76irf1g0xkh6hgpv4qn2swbiyvlazvlrhjfyxb9bckq9"))))
    (version "6.4.99-git")
    (source (local-file "../.." "ogs-checkout"
                      #:recursive? #t
                      #:select? vcs-file?
                      ))
    ; (source (local-file %source-dir #:recursive? #t))
    (build-system cmake-build-system)
    (arguments
     `(#:build-type "Release"
       #:configure-flags (list
                        ;   (string-append "-DOGS_VERSION=" version) ; does not work
                          "-DGUIX_BUILD=ON"
                          "-DOGS_BUILD_TESTING=OFF"
                          "-DOGS_USE_EIGEN_UNSUPPORTED=OFF" ; Eigen 3.4.0
                          "-DOGS_INSTALL_DEPENDENCIES=OFF"  ; handled by guix
                          "-DOGS_CPU_ARCHITECTURE=OFF"      ; enable guix --tune
                          "-DOGS_VERSION=6.4.99"
                          )
       #:cmake ,cmake)) ;for newer CMake version
    (inputs (list boost
                  eigen
                  fmt
                  googletest
                  hdf5
                  json-modern-cxx
                  libxml2
                  pybind11-2.10.4
                  python
                  range-v3
                  spdlog
                  zlib
                  vtk))
    (native-inputs (list git ninja googletest nss-certs)) ; TODO: cpm,
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
        `(cons* "-DOGS_USE_PETSC=ON"
                "-DCMAKE_C_COMPILER=mpicc"
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
    (synopsis "OGS with PETSc and SteadyStateDiffusion only (for faster build testing)")))

(define-public vtk-openmpi
  (package
    (inherit vtk)
    (name "vtk-openmpi")
    (inputs (modify-inputs (package-inputs vtk)
              (prepend hdf5-parallel-openmpi openmpi)))
    (arguments
     (substitute-keyword-arguments (package-arguments vtk)
       ((#:configure-flags flags)
        `(cons* "-DVTK_MODULE_ENABLE_VTK_IOParallelXML=YES"
                "-DVTK_MODULE_ENABLE_VTK_ParallelMPI=YES" "-DVTK_USE_MPI=ON"
                ,flags))))
    (synopsis "VTK with OpenMPI support")))

(define pybind11-2.10.4
      (package
        (inherit pybind11)
        (version "2.10.4")
        (source (origin
                  (method git-fetch)
                  (uri (git-reference
                    (url "https://github.com/pybind/pybind11")
                    (commit (string-append "v" version))))
                  (sha256
                   (base32
                    "0rbcfvl7y472sykzdq3vrkw83kar0lpzhk3wq9yj9cdydl8cpfcz"))))))
;; return package
ogs
