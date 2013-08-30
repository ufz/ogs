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


#include <iostream>
#include <cassert>

#include "logog/include/logog.hpp"

namespace NumLib
{

template <class T_MESH_ELEMENT, class T_SHAPE_FUNC, class T_SHAPE_DATA>
NaturalCoordinatesMapping<T_MESH_ELEMENT,T_SHAPE_FUNC,T_SHAPE_DATA>::NaturalCoordinatesMapping(const MeshElementType &ele)
: _ele(&ele)
{
    this->reset(ele);
};

template <class T_MESH_ELEMENT, class T_SHAPE_FUNC, class T_SHAPE_DATA>
void NaturalCoordinatesMapping<T_MESH_ELEMENT,T_SHAPE_FUNC,T_SHAPE_DATA>::reset(const MeshElementType &ele)
{
    _ele = &ele;
};

template <class T_MESH_ELEMENT, class T_SHAPE_FUNC, class T_SHAPE_DATA>
void NaturalCoordinatesMapping<T_MESH_ELEMENT,T_SHAPE_FUNC,T_SHAPE_DATA>::computeMappingMatrices(const double* natural_pt, ShapeDataType &prop, const ShapeFieldType fields) const
{
    //prepare
    const std::size_t dim = _ele->getDimension();
    const std::size_t nnodes = _ele->getNNodes();
    prop.setZero(fields);

    //shape, dshape/dr
    if ( fields & SHAPE_N )
        T_SHAPE_FUNC::computeShapeFunction(natural_pt, prop.N);

    if ( !(fields & SHAPE_DNDX)) return;

    double* dNdr = prop.dNdr.data();
    T_SHAPE_FUNC::computeGradShapeFunction(natural_pt, dNdr);

    //jacobian: J=[dx/dr dy/dr // dx/ds dy/ds]
    for (std::size_t i_r=0; i_r<dim; i_r++) {
        for (std::size_t j_x=0; j_x<dim; j_x++) {
            for (std::size_t k=0; k<nnodes; k++) {
                prop.J(i_r,j_x) += prop.dNdr(i_r, k) * _ele->getNode(k)->getCoords()[j_x];
            }
        }
    }

    //|J|, J^-1, dshape/dx
    prop.detJ = prop.J.determinant();
    if (prop.detJ>.0) {
        prop.invJ = prop.J.inverse();
        prop.dNdx = prop.invJ * prop.dNdr;
    } else {
        ERR("***error: det_j=%e is not positive.\n", prop.detJ);
    }
};

template <class T_MESH_ELEMENT, class T_SHAPE_FUNC, class T_SHAPE_DATA>
void NaturalCoordinatesMapping<T_MESH_ELEMENT,T_SHAPE_FUNC,T_SHAPE_DATA>::mapToPhysicalCoordinates(const ShapeDataType &prop, double* physical_pt) const
{
    const std::size_t dim = _ele->getDimension();
    const std::size_t nnodes = _ele->getNNodes();

    for (std::size_t i=0; i<dim; i++) {
        physical_pt[i] = .0;
        for (std::size_t j=0; j<nnodes; j++)
            physical_pt[i] += prop.N(j) * _ele->getNode(j)->getCoords()[i];
    }
}

template <class T_MESH_ELEMENT, class T_SHAPE_FUNC, class T_SHAPE_DATA>
void NaturalCoordinatesMapping<T_MESH_ELEMENT,T_SHAPE_FUNC,T_SHAPE_DATA>::mapToNaturalCoordinates(const ShapeDataType &prop, const double* physical_pt, double* natural_pt) const
{
    const std::size_t dim = _ele->getDimension();
    const std::size_t nnodes = _ele->getNNodes();

    // calculate dx which is relative coordinates from the element center
    std::vector<double> dx(dim, .0);
    // x_avg = sum_i {x_i} / n
    for (std::size_t i=0; i<dim; i++)
        for (std::size_t j=0; j<nnodes; j++)
            dx[i] += _ele->getNode(j)->getCoords()[i];
    for (std::size_t i=0; i<dim; i++)
        dx[i] /= (double)nnodes;
    // dx = pt - x_avg
    for (std::size_t i=0; i<dim; i++)
        dx[i] = physical_pt[i] - dx[i];

    // r = invJ^T * dx
    for (std::size_t i=0; i<dim; i++) {
        natural_pt[i] = 0.0;
        for (std::size_t j=0; j<dim; j++)
            natural_pt[i] += prop.invJ(j * dim, i) * dx[j];
    }

}

} //end
