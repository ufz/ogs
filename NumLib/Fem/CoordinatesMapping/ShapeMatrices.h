/**
 * \author Norihiro Watanabe
 * \date   2013-08-13
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <ostream>

#include <Eigen/Core>

namespace NumLib
{

/**
  * \brief Shape matrix type to be calculated
  *
  */
enum class ShapeMatrixType
{
    N,      ///< calculates N
    DNDR,   ///< calculates dNdr
    N_J,    ///< calculates N, dNdr, J, and detJ
    DNDR_J, ///< calculates dNdr, J, and detJ
    DNDX,   ///< calculates dNdr, J, detJ, invJ, and dNdx
    ALL     ///< calculates all
};

/**
 * \brief Coordinates mapping matrices at particular location
 *
 * \tparam T_N      Vector type for shape functions
 * \tparam T_DNDR   Matrix type for gradient of shape functions in natural coordinates
 * \tparam T_J      Jacobian matrix type
 * \tparam T_DNDX   Matrix type for gradient of shape functions in physical coordinates
 */
template <class T_N, class T_DNDR, class T_J, class T_DNDX>
struct ShapeMatrices
{
    using ShapeType = T_N;
    using DrShapeType = T_DNDR;
    using JacobianType = T_J;
    using DxShapeType = T_DNDX;

    ShapeType N;        ///< Vector of shape functions, N(r)
    DrShapeType dNdr;   ///< Matrix of gradient of shape functions in natural coordinates, dN(r)/dr
    JacobianType J;     ///< Jacobian matrix, J=dx/dr
    double detJ;        ///< Determinant of the Jacobian
    JacobianType invJ;  ///< Inverse matrix of the Jacobian
    DxShapeType dNdx;   ///< Matrix of gradient of shape functions in physical coordinates, dN(r)/dx
    double integralMeasure;

    /** Not default constructible, dimensions always must be given.
     *
     * The default constructor has been deleted explicitly, because with
     * dynamically allocated matrices it is rather easy to forget the
     * required <tt>resize()</tt> call. Note: the <tt>resize()</tt> member
     * is also deleted now.
     */
    ShapeMatrices() = delete;

    /**
     * Initialize matrices and vectors
     *
     * @param local_dim  Spatial dimension of the element e.g. 1 for line, 2 for
     *                   quad, and 3 hex etc.
     * @param global_dim Spatial dimension of the element's exterior space e.g.
     *                   3 for a quad representing a boundary of a hex element.
     * @param n_nodes    The number of element nodes
     */
    ShapeMatrices(std::size_t local_dim, std::size_t global_dim,
                  std::size_t n_nodes)
        : N(n_nodes),
          dNdr(local_dim, n_nodes),
          J(local_dim, local_dim),
          detJ(.0),
          invJ(local_dim, local_dim),
          dNdx(global_dim, n_nodes),
          integralMeasure(0.0)
    {
        setZero();
    }

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

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
}; // ShapeMatrices


} // NumLib

#include "ShapeMatrices-impl.h"
