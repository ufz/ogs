/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef TEST_MATRIXTOOLS_H_
#define TEST_MATRIXTOOLS_H_

namespace FeTestData
{

// copy matrix entries in upper triangle to lower triangle
template <class T_MATRIX, typename ID_TYPE=signed>
inline void copyUpperToLower(const ID_TYPE dim, T_MATRIX &m)
{
    for (ID_TYPE i=0; i<dim; i++)
        for (ID_TYPE j=0; j<i; j++)
            m(i,j) = m(j,i);
}

// set an identity matrix
template <class T_MATRIX, typename ID_TYPE=signed>
inline void setIdentityMatrix(unsigned dim, T_MATRIX &m)
{
    for (unsigned i=0; i<dim; i++)
        for (unsigned j=0; j<dim; j++)
            m(i,j) = 0.0;
    for (unsigned i=0; i<dim; i++)
        m(i,i) = 1.0;
}

} // namespace

#endif

