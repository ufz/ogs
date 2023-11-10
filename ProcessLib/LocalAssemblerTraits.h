/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cassert>
#include <type_traits>

#include "NumLib/Fem/ShapeMatrixPolicy.h"

namespace ProcessLib
{

namespace detail
{

/*! This class makes it possible to handle Eigen's block expressions
 * transparently, both for fixed-size and dynamically allocated Eigen matrices.
 *
 * \tparam ShpPol   the shape matrix policy used
 * \tparam NNodes   the number of nodes for the FEM element
 * \tparam NodalDOF the number of d.o.f. per node.
 *                  If zero, this struct's methods will determine the size of
 *                  the generated Eigen::block's at runtime. If non-zero, the
 *                  size will be determined at compile-time.
 * \tparam Dim      global spatial dimension
 */
template <typename ShpPol, unsigned NNodes, unsigned NodalDOF, int Dim>
struct LocalAssemblerTraitsFixed
{
private:
    template <int N, int M>
    using Matrix = typename ShpPol::template MatrixType<N, M>;
    template <int N>
    using Vector = typename ShpPol::template VectorType<N>;

public:
    using ShapeMatrices = typename ShpPol::ShapeMatrices;

    //! Square matrix of the given space dimensionality
    using MatrixDimDim = Matrix<Dim, Dim>;

    // TODO That only works if all process variables are single component
    //      and use the same shape functions.
    //! Local matrix for the given number of d.o.f.\ per node and number of
    //! integration points
    using LocalMatrix = Matrix<NNodes * NodalDOF, NNodes * NodalDOF>;
    //! Local vector for the given number of d.o.f.\ per node and number of
    //! integration points
    using LocalVector = Vector<NNodes * NodalDOF>;

    //! Local vector for one component of one process variable.
    //! The size is the number of nodes in the element.
    using Vector1Comp = Vector<NNodes>;

    //! Laplace matrix for the given space dimensionality and number of d.o.f
    //! per node.
    using LaplaceMatrix = Matrix<Dim * NodalDOF, Dim * NodalDOF>;

    //! Get a block \c Dim x \c Dim whose upper left corner is at \c top and \c
    //! left.
    template <typename Mat>
    static auto blockDimDim(Mat const& mat, unsigned top, unsigned left,
                            unsigned nrows, unsigned ncols)
    {
        if constexpr (NodalDOF != 0)
        {
            assert(nrows == Dim && ncols == Dim);
            (void)nrows;
            (void)ncols;
            return mat.template block<Dim, Dim>(top, left);
        }
        else
        {
            assert(nrows == ncols);
            return mat.block(top, left, nrows, ncols);
        }
    }

    //! Get a block \c Dim x \c Dim whose upper left corner is at \c top and \c
    //! left.
    template <typename Mat>
    static auto blockDimDim(Mat& mat, unsigned top, unsigned left,
                            unsigned nrows, unsigned ncols)
    {
        if constexpr (NodalDOF != 0)
        {
            assert(nrows == Dim && ncols == Dim);
            (void)nrows;
            (void)ncols;
            return mat.template block<Dim, Dim>(top, left);
        }
        else
        {
            assert(nrows == ncols);
            return mat.block(top, left, nrows, ncols);
        }
    }

    //! Get a block \c NNodes x \c NNodes whose upper left corner is at \c top
    //! and \c left.
    template <typename Mat>
    static auto blockShpShp(Mat const& mat, unsigned top, unsigned left,
                            unsigned nrows, unsigned ncols)
    {
        if constexpr (NodalDOF != 0)
        {
            assert(nrows == NNodes && ncols == NNodes);
            (void)nrows;
            (void)ncols;
            return mat.template block<NNodes, NNodes>(top, left);
        }
        else
        {
            assert(nrows == ncols);
            return mat.block(top, left, nrows, ncols);
        }
    }

    //! Get a block \c NNodes x \c NNodes whose upper left corner is at \c top
    //! and \c left.
    template <typename Mat>
    static auto blockShpShp(Mat& mat, unsigned top, unsigned left,
                            unsigned nrows, unsigned ncols)
    {
        if constexpr (NodalDOF != 0)
        {
            assert(nrows == NNodes && ncols == NNodes);
            (void)nrows;
            (void)ncols;
            return mat.template block<NNodes, NNodes>(top, left);
        }
        else
        {
            assert(nrows == ncols);
            return mat.block(top, left, nrows, ncols);
        }
    }

    //! Get a block \c NNodes x 1 starting at the \c top'th row.
    template <typename Vec>
    static auto blockShp(Vec const& vec, unsigned top, unsigned nrows)
    {
        if constexpr (NodalDOF != 0)
        {
            assert(nrows == NNodes);
            (void)nrows;
            return vec.template block<NNodes, 1>(top, 0);
        }
        else
        {
            return vec.block(top, 0, nrows, 1);
        }
    }

    //! Get a block \c NNodes x 1 starting at the \c top'th row.
    template <typename Vec>
    static auto blockShp(Vec& vec, unsigned top, unsigned nrows)
    {
        if constexpr (NodalDOF != 0)
        {
            assert(nrows == NNodes);
            (void)nrows;
            return vec.template block<NNodes, 1>(top, 0);
        }
        else
        {
            return vec.block(top, 0, nrows, 1);
        }
    }
};

}  // namespace detail

#ifndef OGS_EIGEN_DYNAMIC_SHAPE_MATRICES

template <typename ShpPol, unsigned NNodes, unsigned NodalDOF, int Dim>
using LocalAssemblerTraits =
    detail::LocalAssemblerTraitsFixed<ShpPol, NNodes, NodalDOF, Dim>;

static_assert(OGS_EIGEN_DYNAMIC_SHAPE_MATRICES_FLAG == 0,
              "Inconsistent use of the macro OGS_EIGEN_DYNAMIC_SHAPE_MATRICES."
              " Maybe you forgot to include some header file.");
#else

template <typename ShpPol, unsigned NNodes, unsigned NodalDOF, int Dim>
using LocalAssemblerTraits =
    detail::LocalAssemblerTraitsFixed<ShapeMatrixPolicyType<void, 0>, 0, 0, 0>;

static_assert(OGS_EIGEN_DYNAMIC_SHAPE_MATRICES_FLAG == 1,
              "Inconsistent use of the macro OGS_EIGEN_DYNAMIC_SHAPE_MATRICES."
              " Maybe you forgot to include some header file.");
#endif

}  // namespace ProcessLib
