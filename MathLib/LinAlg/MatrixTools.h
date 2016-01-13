/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATHLIB_MATRIXTOOLS_H_
#define MATHLIB_MATRIXTOOLS_H_

#include <type_traits>

#ifdef OGS_USE_EIGEN
#include <Eigen/Eigen>
#endif


namespace MathLib
{

#ifdef OGS_USE_EIGEN
namespace details
{

template <class T, int N, typename = void>
struct Determinant
{
    static_assert(N != Eigen::Dynamic, "MathLib::Determinant() for a fixed-size matrix  was called for a dynamic one");
    static double eval(T const& mat)
    {
        return mat.determinant();
    }
};

template <class T, int N>
struct Determinant<T, N, typename std::enable_if<N == Eigen::Dynamic>::type>
{
    static_assert(N == Eigen::Dynamic, "MathLib::Determinant() for a dynamic-size matrix  was called for a static one");
    static double eval(T const& mat)
    {
        if (mat.rows()==1)
            return mat(0,0);
        else if (mat.rows()==2)
            return mat(0,0) * mat(1,1) - mat(0,1) * mat(1,0);
        else if (mat.rows()==3)
            return mat(0,0) * (mat(1,1) * mat(2,2) - mat(2,1) * mat(1,2))
                     + mat(2,0) * (mat(0,1) * mat(1,2) - mat(1,1) * mat(0,2))
                     + mat(1,0) * (mat(0,2) * mat(2,1) - mat(2,2) * mat(0,1));
        else
            return mat.determinant();
    }
};


template <class T, int N, typename = void>
struct Inverse
{
    static_assert(N != Eigen::Dynamic, "MathLib::Inverse() for a fixed-size matrix  was called for a dynamic one");
    static void eval(T const& mat, double /*detJ*/, T &result)
    {
        result.noalias() = mat.inverse();
    }
};

template <class T, int N>
struct Inverse<T, N, typename std::enable_if<N == Eigen::Dynamic>::type>
{
    static_assert(N == Eigen::Dynamic, "MathLib::Inverse() for a dynamic-size matrix  was called for a static one");
    static void eval(T const& mat, double detJ, T &result)
    {
        if (mat.rows()==1) {
            result(0,0) = 1./detJ;
        } else if (mat.rows()==2) {
            result(0,0) = mat(1,1);
            result(0,1) = -mat(0,1);
            result(1,0) = -mat(1,0);
            result(1,1) = mat(0,0);
            result *= 1./detJ;
        } else if (mat.rows()==3) {
            result(0,0) =  mat(1,1) * mat(2,2) - mat(2,1) * mat(1,2);
            result(0,1) =  mat(0,2) * mat(2,1) - mat(0,1) * mat(2,2);
            result(0,2) =  mat(0,1) * mat(1,2) - mat(0,2) * mat(1,1);
            result(1,0) =  mat(1,2) * mat(2,0) - mat(2,2) * mat(1,0);
            result(1,1) =  mat(0,0) * mat(2,2) - mat(2,0) * mat(0,2);
            result(1,2) =  mat(0,2) * mat(1,0) - mat(1,2) * mat(0,0);
            result(2,0) =  mat(1,0) * mat(2,1) - mat(2,0) * mat(1,1);
            result(2,1) =  mat(0,1) * mat(2,0) - mat(2,1) * mat(0,0);
            result(2,2) =  mat(0,0) * mat(1,1) - mat(1,0) * mat(0,1);
            result *= 1./detJ;
        } else {
            result.noalias() = mat.inverse();
        }
    }
};

} // details

/// returns determinant of a given Eigen dense matrix
template <typename Derived>
double determinant(Eigen::MatrixBase<Derived> const& mat)
{
    using MatrixType = Eigen::MatrixBase<Derived>;
    return details::Determinant<MatrixType, MatrixType::SizeAtCompileTime>::eval(mat);
}

/// computes inverse of a given Eigen dense matrix
/// \param mat     a dense matrix whose inverse is computed
/// \param det_mat determinant of the given matrix
/// \param result  inverted matrix whose memory should be allocated beforehand
template <typename Derived>
void inverse(Eigen::MatrixBase<Derived> const& mat, double det_mat, Eigen::MatrixBase<Derived> & result)
{
    using MatrixType = Eigen::MatrixBase<Derived>;
    return details::Inverse<MatrixType, MatrixType::SizeAtCompileTime>::eval(mat, det_mat, result);
}


#endif

} // namespace

#endif /* MATHLIB_MATRIXTOOLS_H_ */
