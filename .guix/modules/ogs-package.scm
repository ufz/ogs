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
  #:use-module (guix build-system copy)
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
      (const #t)))            ; not in a Git checkout

(define-public ogs
  (package
    (name "ogs")
    (version "6.4.99-git")
    (source (local-file "../.." "ogs-checkout"
                      #:recursive? #t
                      #:select? vcs-file?
                      ))
    (build-system cmake-build-system)
    (arguments
     `(#:build-type "Release"
       #:configure-flags (list
        ; passing variables works like this
        ;                  ,(string-append "-DOGS_VERSION=" version)
        ; TODO: but it is not overwritten in sub-packages...
                          "-DGUIX_BUILD=ON"
                          "-DOGS_BUILD_TESTING=OFF"
                          "-DOGS_USE_EIGEN_UNSUPPORTED=OFF" ; Eigen 3.4.0
                          "-DOGS_INSTALL_DEPENDENCIES=OFF"  ; handled by guix
                          "-DOGS_CPU_ARCHITECTURE=OFF"      ; enable guix --tune
                          )
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
    (native-inputs (list git ninja nss-certs)) ; TODO: cpm
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

; #### Releases ####
(define-public ogs-6.4.4
  (package
    (inherit ogs)
    (name "ogs")
    (version "6.4.4")
    (source (origin
              (method git-fetch)
              (uri (git-reference
                    (url "https://gitlab.opengeosys.org/ogs/ogs.git")
                    (commit "d4ca7e627f2fc012bfe434649b797e78e5c2a8f1")
                    (recursive? #t)))
              (file-name (git-file-name name version))
              (sha256
               ;; Update with `guix hash -rx .`, make sure to have submodules updated!
               (base32
                "1f6mcbjx76irf1g0xkh6hgpv4qn2swbiyvlazvlrhjfyxb9bckq9"))))
    (synopsis "OGS 6.4.4 release")))

(define-public ogs-petsc-6.4.4
  (package
    (inherit ogs-petsc)
    (name "ogs-petsc")
    (version "6.4.4")
    (source (origin
              (method git-fetch)
              (uri (git-reference
                    (url "https://gitlab.opengeosys.org/ogs/ogs.git")
                    (commit "d4ca7e627f2fc012bfe434649b797e78e5c2a8f1")
                    (recursive? #t)))
              (file-name (git-file-name name version))
              (sha256
               ;; Update with `guix hash -rx .`, make sure to have submodules updated!
               (base32
                "1f6mcbjx76irf1g0xkh6hgpv4qn2swbiyvlazvlrhjfyxb9bckq9"))))
    (synopsis "OGS 6.4.4 with PETSc release")))

; #### Dependencies ####
(define-public vtk-openmpi
  (package
    (inherit vtk)
    (name "vtk-openmpi")
    (inputs (modify-inputs (package-inputs vtk)
              (prepend hdf5-parallel-openmpi openmpi)))
    (arguments
     (substitute-keyword-arguments (package-arguments vtk)
       ((#:configure-flags flags)
        #~(append '("-DVTK_MODULE_ENABLE_VTK_IOParallelXML=YES"
                "-DVTK_MODULE_ENABLE_VTK_ParallelMPI=YES" "-DVTK_USE_MPI=ON")
                #$flags))))
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

(define tetgen
      (package
        (name "tetgen")
        (synopsis "A Quality Tetrahedral Mesh Generator and a 3D Delaunay Triangulator")
        (license license:agpl3)
        (description "TetGen is a program to generate tetrahedral meshes of any 3D polyhedral domains. TetGen generates exact constrained Delaunay tetrahedralizations, boundary conforming Delaunay meshes, and Voronoi partitions.")
        (home-page "http://www.tetgen.org/")
        (version "1.5.1-1")
        (source (origin
                  (method git-fetch)
                  (uri (git-reference
                    (url "https://github.com/ufz/tetgen")
                    (commit version)))
                  (sha256
                   (base32
                    "1xp1qibm0q4z5qx0h178qpas3n7pqbladkxdalq9j4l98hdws46j"))))
        (build-system cmake-build-system)
        (arguments
            `(#:tests? #f)
        )))

(define tclap
      (package
        (name "tclap")
        (synopsis "Templatized Command Line Argument Parser")
        (license license:expat)
        (description "This is a simple C++ library that facilitates parsing command line arguments in a type independent manner.")
        (home-page "https://sourceforge.net/p/tclap/discussion/")
        (version "1.2.4-1")
        (source (origin
                  (method git-fetch)
                  (uri (git-reference
                    (url "https://github.com/ufz/tclap")
                    (commit version)))
                  (sha256
                   (base32
                    "0bijzfc9c8zny3m74y53i8m3f41kd8klcnmr9chy536syr9vdr5p"))))
        (build-system cmake-build-system)
        (arguments
            `(#:tests? #f)
        )))

(define-public iphreeqc
      (package
        (name "iphreeqc")
        (synopsis "Modules Based on the Geochemical Model PHREEQC for use in scripting and programming languages")
        (license license:public-domain)
        (description "")
        (home-page "https://www.usgs.gov/software/phreeqc-version-3")
        (version "3.5.0-3")
        (source (origin
                  (method git-fetch)
                  (uri (git-reference
                    (url "https://github.com/ufz/iphreeqc")
                    (commit version)))
                  (sha256
                   (base32
                    "12wiqyzpzx89k9c7q07w4ypnppvi6s88k6jjsnlnvaxfafyvrbw3"))))
        (build-system cmake-build-system)
        (arguments
            `(#:tests? #f)
        )))

(define xmlpatch
      (package
        (name "xmlpatch")
        (synopsis "An XML Patch library")
        (license license:lgpl2.1)
        (description "")
        (home-page "https://xmlpatch.sourceforge.net")
        (version "0.4.3")
        (source (origin
                  (method git-fetch)
                  (uri (git-reference
                    (url "https://gitlab.opengeosys.org/ogs/libs/xmlpatch")
                    (commit (string-append "v" version))))
                  (sha256
                   (base32
                    "0872g9w1jd5r4c5a1s8ga4x1plg608b7rxyqjs6zv8ghjq9qlkvg"))))
        (build-system cmake-build-system)
        (inputs (list libxml2))
        (arguments
            `(#:tests? #f)
        )))

(define exprtk
      (package
        (name "exprtk")
        (home-page "https://www.partow.net/programming/exprtk/index.html")
        (synopsis "C++ Mathematical Expression Parsing And Evaluation Library")
        (description "")
        (license license:expat)
        (version "0.0.2")
        (source (origin
                  (method git-fetch)
                  (uri (git-reference
                    (url "https://github.com/ArashPartow/exprtk.git")
                    (commit version)))
                  (sha256
                   (base32
                    "1w92qlfjpcan38d88fak3avq81lkcpai5mqpbvrsfv04mi5nfpk5"))))
        (build-system copy-build-system)
        (arguments
         '(#:install-plan '(("exprtk.hpp" "include/")))
        )))

;; return package
ogs

;; TODO: add this to web page
;  ## Using the ogs repo as a Guix channel (as a user of ogs)
;
;  Add the following to `~/.config/guix/channels.scm`:
;
;  ```scheme
;  (append (list (channel
;                  (name 'ogs)
;                  (url "https://gitlab.opengeosys.org/ogs/ogs.git")
;                  (branch "master"))) %default-channels)
;  ```
;
;  Run `guix pull`.
;
;  Then you can install the ogs package:
;
;  ```bash
;  guix install ogs@6.4.4
;  # OR, e.g.
;  guix install ogs-petsc@6.4.4
;  ```
