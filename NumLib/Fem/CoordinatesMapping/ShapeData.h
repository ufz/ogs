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


#ifndef SHAPEDATA_H_
#define SHAPEDATA_H_

namespace NumLib
{

/**
 * \brief Coordinate mapping matrices at particular location
 *
 * Mapping is done between physical coordinates (x,y,z) and natural coordinates (r,s,t)
 */
template <class T_N, class T_DN, class T_J>
struct ShapeData
{
    typedef T_N ShapeType;
    typedef T_DN DShapeType;
    typedef T_J JacobianType;

    /// shape function N(r)
    ShapeType N;

    /// gradient of shape functions, dN(r)/dr
    DShapeType dNdr;

    /// gradient of shape functions, dN(r)/dx
    DShapeType dNdx;

    /// Jacobian matrix, J=dx/dr
    JacobianType J;

    /// inverse of the Jacobian
    JacobianType invJ;

    /// determinant of the Jacobian
    double detJ;

    /**
     *
     * @param dim       Dimension of the physical coordinates
     * @param n_nodes   The number of element nodes
     */
    ShapeData(std::size_t dim, std::size_t n_nodes)
    : N(n_nodes), dNdr(dim, n_nodes),
      dNdx(dim, n_nodes), J(dim, dim),
      invJ(dim, dim), detJ(.0)
    {}

    ~ShapeData() {}

    void setZero()
    {
        setZero(N);
        setZero(dNdr);
        setZero(dNdx);
        setZero(J);
        setZero(invJ);
        detJ = .0;
    }

private:
    template<class T>
    void setZero(T &mat)
    {
        mat.setZero(mat.rows(), mat.cols());
    }

    void setZero(ShapeType &vec)
    {
        vec.setZero(vec.size());
    }

};

}

#endif //SHAPEDATA_H_
