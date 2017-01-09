/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef SHAPEHEX8_H_
#define SHAPEHEX8_H_

#include "MeshLib/Elements/Hex.h"

namespace NumLib
{

/**
 *  Shape function for a 8-nodes hex element in natural coordinates
 *
 * \verbatim
 *       (-1, 1, 1) 7-----------6 ( 1, 1, 1)
 *                 /:          /|
 *                / :         / |
 *               /  :        /  |
 *              /   :       /   |
 *             /    :      /    |
 * (-1,-1, 1) 4-----------5 ( 1,-1, 1)
 *            |     :     |     |
 *       (-1, 1,-1) 3.....|.....2 ( 1, 1,-1)
 *            |    .      |    /
 *            |   .       |   /
 *            |  .        |  /
 *            | .         | /
 *            |.          |/
 * (-1,-1,-1) 0-----------1 ( 1,-1,-1)
 *          0
 * \endverbatim
 */
class ShapeHex8
{
public:
    /**
     * Evaluate the shape function at the given point
     *
     * @param [in]  r   natural coordinates (r,s,t)
     * @param [out] N   a vector of calculated shape functions
     */
    template <class T_X, class T_N>
    static void computeShapeFunction(const T_X &r, T_N &N);

    /**
     * Evaluate derivatives of the shape function at the given point
     *
     * @param [in]  r   natural coordinates (r,s,t)
     * @param [out] dN  a matrix of the derivatives
     */
    template <class T_X, class T_N>
    static void computeGradShapeFunction(const T_X &r, T_N &dN);

    using MeshElement = MeshLib::Hex;
    static const unsigned DIM = MeshElement::dimension;
    static const unsigned NPOINTS = MeshElement::n_all_nodes;
};

}

#include "ShapeHex8-impl.h"

#endif //SHAPEHEX8_H_
