Enables scaling algorithm to equilibrate rows and column norms in matrices
for a LU/ILU factorization if the unsupported Eigen modules are available
through OGS_USE_EIGEN_UNSUPPORTED cmake option.

See following links for detail description:
 - "Note" in https://eigen.tuxfamily.org/dox/classEigen_1_1SparseLU.html
 - "Detailed Description" in https://eigen.tuxfamily.org/dox/unsupported/namespaceEigen.html
 - References : D. Ruiz and B. Ucar, A Symmetry Preserving Algorithm for Matrix Scaling, INRIA Research report RR-7552
