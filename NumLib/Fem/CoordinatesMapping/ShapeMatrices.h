/**
 * \author Norihiro Watanabe
 * \date   2013-08-13
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef SHAPEMATRICES_H_
#define SHAPEMATRICES_H_

#include <iostream>

namespace NumLib
{

/**
  * \brief Shape matrix type to be calculated
  *
  */
enum class ShapeMatrixType
{
    N,      //< calculates N
    DNDR,   //< calculates dNdr
    N_J,    //< calculates N, dNdr, J, and detJ
    DNDR_J, //< caluclates dNdr, J, and detJ
    DNDX,   //< calculates dNdr, J, detJ, invJ, and dNdx
    ALL     //< calculates all
};

/**
 * \brief Coordinates mapping matrices at particular location
 *
 * \tparam T_N      Vector type for shape functions
 * \tparam T_DN     Matrix type for gradient of shape functions
 * \tparam T_J      Jacobian matrix type
 */
template <class T_N, class T_DN, class T_J>
struct ShapeMatrices
{
    typedef T_N ShapeType;
    typedef T_DN DShapeType;
    typedef T_J JacobianType;

    ShapeType N;        ///< Vector of shape functions, N(r)
    DShapeType dNdr;    ///< Matrix of gradient of shape functions in natural coordinates, dN(r)/dr
    JacobianType J;     ///< Jacobian matrix, J=dx/dr
    double detJ;        ///< Determinant of the Jacobian
    JacobianType invJ;  ///< Inverse matrix of the Jacobian
    DShapeType dNdx;    ///< Matrix of gradient of shape functions in physical coordinates, dN(r)/dx

    /** The default constructor is used by fixed-size (at compile-time)
     * matrix/vector types where no resizing of the matrices/vectors is required
     * in this constructor.
     */
    ShapeMatrices() = default;

    /**
     * Initialize matrices and vectors
     *
     * @param dim       Spatial dimension
     * @param n_nodes   The number of element nodes
     */
    ShapeMatrices(std::size_t dim, std::size_t n_nodes)
    : N(n_nodes), dNdr(dim, n_nodes), J(dim, dim), detJ(.0),
      invJ(dim, dim), dNdx(dim, n_nodes)
    {
        this->setZero();
    }

    ~ShapeMatrices() {}

    /// Reinitialize the matrices and vectors. When using dynamic size matrices
    /// and the default constructor of this class memory allocation must be
    /// performed calling this function.
    /// For fixed size matrices no memory reallocation happens and the
    /// matrix/vector sizes must be the same as at construction (given by the
    /// template parameters).
    void resize(std::size_t const dim, std::size_t n_nodes);

    /// reset all data with zero
    void setZero();

    /**
     * reset specified data with zero
     *
     * @tparam T_SHAPE_MATRIX_TYPE    shape matrix types to be initialized
     */
    template <ShapeMatrixType T_SHAPE_MATRIX_TYPE>
    void setZero();

    /**
     * writes the matrix entries into the output stream
     * @param out the output stream
     */
    void write (std::ostream& out) const;
}; // ShapeMatrices


} // NumLib

#include "ShapeMatrices-impl.h"

#endif //SHAPEMATRICES_H_
