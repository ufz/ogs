/**
 * \file
 * \author Thomas Fischer
 * \date   2011-05-06
 * \brief  Definition of the GaussAlgorithm class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef GAUSSALGORITHM_H_
#define GAUSSALGORITHM_H_

#include <vector>


#include "BaseLib/ConfigTree.h"
#include "../Dense/DenseMatrix.h"
#include "TriangularSolve.h"

namespace MathLib {

/**
 * This is a class for the direct solution of (dense) systems of
 * linear equations, \f$A x = b\f$. During the construction of
 * the object the matrix A is factorized in matrices L and U using
 * Gauss-Elimination with partial pivoting (rows are exchanged). In doing so
 * the entries of A change! The solution for a specific
 * right hand side is computed by the method execute().
 */
template <typename MAT_T, typename VEC_T = typename MAT_T::FP_T*>
class GaussAlgorithm
{
public:
    typedef typename MAT_T::FP_T FP_T;
    typedef typename MAT_T::IDX_T IDX_T;

public:
    /**
     * A direct solver for the (dense) linear system \f$A x = b\f$.
     * @param solver_name A name used as a prefix for command line options
     *                    if there are such options available.
     * @param option For some solvers the user can give parameters to the
     * algorithm. GaussAlgorithm has to fulfill the common interface
     * of all solvers of systems of linear equations. For this reason the
     * second argument was introduced.
     */
    GaussAlgorithm(const std::string solver_name = "",
                   BaseLib::ConfigTree const*const option = nullptr)
    {
        (void) solver_name; (void) option; // silence both compiler and doxygen warnings.
    }

    /**
     * Method solves the linear system \f$A x = b\f$ (based on the LU factorization)
     * using forward solve and backward solve.
     * @param A the coefficient matrix
     * @param b at the beginning the right hand side, at the end the solution
     * @param decompose Flag that signals if the LU decomposition should be
     *        performed or not. If the matrix \f$A\f$ does not change, the LU
     *        decomposition needs to be performed once only!
     * @attention The entries of the given matrix will be changed!
     */
    template <typename V>
    void solve (MAT_T& A, V & b, bool decompose = true);

    void solve(MAT_T& A, FP_T* & b, bool decompose = true);

    /**
     * Method solves the linear system \f$A x = b\f$ (based on the LU factorization)
     * using forward solve and backward solve.
     * @param A (input) the coefficient matrix
     * @param b (input) the right hand side
     * @param x (output) the solution
     * @param decompose see documentation of the other solve methods.
     * @attention The entries of the given matrix will be changed!
     */
    void solve(MAT_T& A, VEC_T const& b, VEC_T & x, bool decompose = true);

private:
    // void solve (MAT_T& A, VEC_T const& b, bool decompose);

    void performLU(MAT_T& A);
    /**
     * permute the right hand side vector according to the
     * row permutations of the LU factorization
     * @param b the entries of the vector b are permuted
     */
    template <typename V> void permuteRHS(V & b) const;
    void permuteRHS (VEC_T& b) const;

    //! the permutation of the rows
    std::vector<IDX_T> _perm;
};

} // end namespace MathLib

#include "GaussAlgorithm-impl.h"

#endif /* GAUSSALGORITHM_H_ */
