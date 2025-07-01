(use-modules (guix transformations)
             (guix packages)
             (guix build-system copy)
             (guix build utils)
             (guix gexp)
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
                     (specification->package "autocheck")
                     ;; OGS sources inside container
                     ;; comment out for smaller dev containert without sources
                     ogs-source-package
                    ))
