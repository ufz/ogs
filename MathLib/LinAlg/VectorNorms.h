/**
 * \file
 * \author Thomas Fischer
 * \author Wenqing Wang
 *
 * \date   2011-06-06 -- 2013-12-10
 *
 * \brief  Definition of vector norm functions.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef VECTORNORMS_H_
#define VECTORNORMS_H_

#include <cmath>
#include <string>

#include "MathTools.h"
#include "LinAlgEnums.h"

namespace MathLib
{

// ----------------------------------------------------------------------------
// General norm functions
// ----------------------------------------------------------------------------

/// return L1 norm \f$ \|x\|_1 = \sum_{i=1}^n |x_i| \f$
/// \tparam T   Data type
/// \param  v   a vector object
/// \param  n   size of the vector
template<class T>
inline T norm_1(T const * const v, std::size_t n)
{
    T s = .0;
    for (std::size_t i=0; i<n; i++)
        s += std::abs(v[i]);
    return s;
}

/// return L2 norm (Euclidean norm) \f$ \|x\|_2 = \sqrt{\sum_{i=1}^n x_i^2} \f$
/// \tparam T   Data type
/// \param  v   a vector object
/// \param  n   size of the vector
template<class T>
inline T norm_2(T const * const v, std::size_t n)
{
    return std::sqrt (scalarProduct (v, v, n));
}

/// return maximum norm \f$ \|x\|_\infty = \max(|x_1|, ..., |x_n|) \f$
/// \tparam T   Data type
/// \param  v   a vector object
/// \param  n   size of the vector
template<class T>
inline T norm_max(T const * const v, std::size_t n)
{
    T u_max = .0;
    for (std::size_t i=0; i<n; i++)
        u_max = std::max(u_max, std::abs(v[i]));
    return u_max;
}


// ----------------------------------------------------------------------------
// Norm functions with a single argument
// ----------------------------------------------------------------------------

/// return L1 norm \f$ \|x\|_1 = \sum_{i=1}^n |x_i| \f$
/// \tparam T   Vector type
/// \param  v   a vector object
template<class T_VEC>
inline double norm_1(const T_VEC &v) { return MathLib::norm_1(&v[0], v.size()); }

/// return L2 norm (Euclidean norm) \f$ \|x\|_2 = \sqrt{\sum_{i=1}^n x_i^2} \f$
/// \tparam T   Vector type
/// \param  v   a vector object
template<class T_VEC>
inline double norm_2(const T_VEC &v) { return MathLib::norm_2(&v[0], v.size()); }

/// return maximum norm \f$ \|x\|_\infty = \max(|x_1|, ..., |x_n|) \f$
/// \tparam T   Vector type
/// \param  v   a vector object
template<class T_VEC>
inline double norm_max(const T_VEC &v) { return MathLib::norm_max(&v[0], v.size()); }

/// return absolute value
/// \param  v   value
template<>
inline double norm_1(const double &v) { return std::abs(v); }

/// return absolute value
/// \param  v   a vector object
template<>
inline double norm_2(const double &v) { return norm_1(v); }

/// return absolute value
/// \param  v   a vector object
template<>
inline double norm_max(const double &v) { return norm_1(v); }

/**
 * return norm
 *
 * \tparam T            vector type
 * \tparam v            a vector object
 * \tparam normType     norm type (default 2-norm)
 * \return calculated norm. If the give norm type is invalid, a negative value is returned.
 */
template<class T_VEC>
inline double norm(const T_VEC &v, VecNormType normType = VecNormType::NORM2)
{
    double norm = -1;
    switch (normType)
    {
    case VecNormType::NORM1:
        norm = norm_1(v);
        break;
    case VecNormType::NORM2:
        norm = norm_2(v);
        break;
    case VecNormType::INFINITY_N:
        norm = norm_max(v);
        break;
    default:
        break;
    }
    return norm;
}


} // end namespace MathLib

#endif /* VECTORNORMS_H_ */
