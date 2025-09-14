(use-modules (guix transformations)
             (guix packages)
             (guix build-system copy)
             (guix build utils)
             (guix gexp)
             (guix git-download)
             (guix utils)
             (gnu packages base)
             ((guix licenses)
              #:prefix license:))

(define current-dir
  (getcwd))

(define current-source
  (local-file current-dir "my-source"
              #:recursive? #t))

(define ogs-source-package
  (package
    (name "ogs-source")
    (version "1.0")
    (source
     current-source)
    (build-system copy-build-system)
    (synopsis "Current ogs source directory")
    (description "Current ogs source directory for development")
    (home-page "http://www.opengeosys.org")
    (license license:bsd-3)))

(define transform1
  (options->transformation `((with-commit . "eigen=9000b3767770f6dd0f4cfb12f4e19c71921885a4")
                             (without-tests . "eigen")
                             (with-configure-flag . "vtk=-DVTK_MODULE_USE_EXTERNAL_VTK_eigen=OFF")
                             (with-source . "vtk=https://www.vtk.org/files/release/9.3/VTK-9.3.1.tar.gz"))))

(define autocheck
  (package
    (name "autocheck")
    (version "0.0.1")
    (source
     (origin
       (method git-fetch)
       (uri (git-reference
             (url "https://github.com/ufz/autocheck")
             (commit "e388ecbb31c49fc2724c8d0436da313b6edca7fd")))
       (file-name (git-file-name name version))
       (sha256
        (base32 "0skw2mk95d8mizzhfm29d5ka71a3y9c772c830yiwhanfpyslql0"))))
    (build-system copy-build-system)
    (arguments
     '(#:install-plan '(("include/autocheck" "include/"))))
    (synopsis "QuickCheck and SmallCheck clones for C++")
    (description "QuickCheck and SmallCheck clones for C++")
    (home-page "https://github.com/thejohnfreeman/autocheck")
    (license license:isc)))

(packages->manifest (list
                     ;; Base packages
                     (specification->package "bash")
                     (specification->package "glibc-locales")
                     (specification->package "nss-certs")

                     ;; Common command line tools
                     (specification->package "coreutils")
                     (specification->package "grep")
                     (specification->package "which")
                     (specification->package "wget")
                     (specification->package "sed")

                     ;; Toolchain
                     (specification->package "gcc-toolchain@13.3.0")
                     (specification->package "ninja")
                     (specification->package "cmake")
                     (specification->package "git")
                     (specification->package "pkg-config")

                     ;; Dependencies
                     (transform1 (specification->package "vtk"))
                     (transform1 (specification->package "vtkdiff"))
                     (transform1 (specification->package "eigen"))
                     (specification->package "boost")
                     (specification->package "exprtk")
                     (specification->package "hdf5")
                     (specification->package "iphreeqc")
                     (specification->package "nlohmann-json")
                     (specification->package "libxml2")
                     (specification->package "mgis")
                     (specification->package "netcdf")
                     (specification->package "netcdf-cxx4")
                     (specification->package "pybind11@2.10.4")
                     (specification->package "python")
                     (specification->package "range-v3")
                     (specification->package "spdlog")
                     (specification->package "tclap")
                     (specification->package "tetgen")
                     (specification->package "tfel")
                     (specification->package "zlib")
                     (specification->package "xmlpatch")
                     (specification->package "metis")
                     (specification->package "googletest")
                     autocheck
                     ;; OGS sources inside container
                     ;; comment out for smaller dev containert without sources
                     ogs-source-package))
