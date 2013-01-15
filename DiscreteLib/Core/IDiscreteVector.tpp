/**
 * \file   IDiscreteVector.tpp
 * \author Norihiro Watanabe
 * \date   2012-08-03
 * \brief  Helper macros.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef IDISCRETEVECTOR_TPP_
#define IDISCRETEVECTOR_TPP_

#include <vector>

namespace DiscreteLib
{

template <typename T>
typename IDiscreteVector<T>::MyVectorType& IDiscreteVector<T>::operator= (const MyVectorType &src)
{
    for (std::size_t i=getRangeBegin(); i<getRangeEnd(); i++)
        (*this)[i] = src[i];
    return *this;
}

template <typename T>
void IDiscreteVector<T>::operator+= (const MyVectorType& v)
{
    for (std::size_t i=getRangeBegin(); i<getRangeEnd(); i++)
        (*this)[i] += v[i];
}
template <typename T>
void IDiscreteVector<T>::operator-= (const MyVectorType& v)
{
    for (std::size_t i=getRangeBegin(); i<getRangeEnd(); i++)
        (*this)[i] -= v[i];
}

template <typename T>
typename IDiscreteVector<T>::MyVectorType& IDiscreteVector<T>::operator= (T v)
{
    for (std::size_t i=getRangeBegin(); i<getRangeEnd(); i++)
        (*this)[i] = v;
    return *this;
}

template <typename T>
void IDiscreteVector<T>::addSubvector(const std::vector<std::size_t> &pos, T* local_v)
{
    for (std::size_t i=0; i<pos.size(); ++i) {
        if (pos[i]==index_npos) continue;
        (*this)[pos[i]] += local_v[i];
    }
}

template <typename T>
void IDiscreteVector<T>::setSubvector(const std::vector<std::size_t> &pos, T v)
{
    for (std::size_t i=0; i<pos.size(); ++i) {
        if (pos[i]==index_npos) continue;
        (*this)[pos[i]] = v;
    }
}

} // end

#endif //IDISCRETEVECTOR_TPP_
