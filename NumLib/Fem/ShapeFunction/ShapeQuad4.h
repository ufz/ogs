/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-08-13
 * \brief
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef SHAPEQUAD4_H_
#define SHAPEQUAD4_H_

namespace NumLib
{

/**
 *  Shape function for a quadrilateral element of four nodes in natural coordinates
 *
 *  2 (-1,1)     1 (1,1)
 *     *--------*
 *     |        |
 *     |        |
 *     |        |
 *     *--------*
 *  3 (-1,-1)    4 (1,-1)
 */
class ShapeQuad4
{
public:
    /**
     * Evaluate the shape function at the given point
     *
     * @param [in]  r    point coordinates
     * @param [out] N   a vector of calculated shape function.
     */
    template <class T_X, class T_N>
    static void computeShapeFunction(const T_X &r, T_N &N4);

    /**
     * Evaluate derivatives of the shape function at the given point
     *
     * @param [in]  r    point coordinates
     * @param [out] dN  a matrix of the derivatives
     */
    template <class T_X, class T_N>
    static void computeGradShapeFunction(const T_X &r, T_N &dN8);
};

}

#include "ShapeQuad4.tpp"

#endif //SHAPEQUAD4_H_
