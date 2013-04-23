/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Implementation of the VectorBase class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef VECTORBASE_TPP_
#define VECTORBASE_TPP_

namespace MathLib
{

template <typename T>
VectorBase<T>& VectorBase<T>::operator= (const VectorBase<T> &src)
{
    for (unsigned i=getRangeBegin(); i<getRangeEnd(); i++)
        this->set(i, src[i]);
    return *this;
}

template <typename T>
void VectorBase<T>::operator+= (const VectorBase<T>& v)
{
    for (unsigned i=getRangeBegin(); i<getRangeEnd(); i++)
        this->add(i, v[i]);
}

template <typename T>
void VectorBase<T>::operator-= (const VectorBase<T>& v)
{
    for (unsigned i=getRangeBegin(); i<getRangeEnd(); i++)
        this->add(i, -v[i]);
}

template <typename T>
VectorBase<T>& VectorBase<T>::operator= (T val)
{
    for (unsigned i=getRangeBegin(); i<getRangeEnd(); i++)
        this->set(i, val);
    return *this;
}

} // MathLib

#endif /* VECTORBASE_TPP_ */
