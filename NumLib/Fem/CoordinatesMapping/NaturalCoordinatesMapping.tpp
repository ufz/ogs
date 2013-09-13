/**
 * \author Norihiro Watanabe
 * \date   2013-08-13
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

namespace detail
{

template <ShapeMatrixType FIELD_TYPE> struct FieldType {};

template <class T_MESH_ELEMENT, class T_SHAPE_FUNC, class T_SHAPE_MATRICES>
inline void computeMappingMatrices(
        const T_MESH_ELEMENT &/*ele*/,
        const double* natural_pt,
        T_SHAPE_MATRICES &shapemat,
        FieldType<ShapeMatrixType::N>)
{
    T_SHAPE_FUNC::computeShapeFunction(natural_pt, shapemat.N);
};

template <class T_MESH_ELEMENT, class T_SHAPE_FUNC, class T_SHAPE_MATRICES>
inline void computeMappingMatrices(
        const T_MESH_ELEMENT &/*ele*/,
        const double* natural_pt,
        T_SHAPE_MATRICES &shapemat,
        FieldType<ShapeMatrixType::DNDR>)
{
    double* const dNdr = shapemat.dNdr.data();
    T_SHAPE_FUNC::computeGradShapeFunction(natural_pt, dNdr);
};

template <class T_MESH_ELEMENT, class T_SHAPE_FUNC, class T_SHAPE_MATRICES>
inline void computeMappingMatrices(
        const T_MESH_ELEMENT &ele,
        const double* natural_pt,
        T_SHAPE_MATRICES &shapemat,
        FieldType<ShapeMatrixType::DNDR_J>)
{
    computeMappingMatrices<T_MESH_ELEMENT, T_SHAPE_FUNC, T_SHAPE_MATRICES>
        (ele, natural_pt, shapemat, FieldType<ShapeMatrixType::DNDR>());

    const std::size_t dim = T_MESH_ELEMENT::dimension;
    const std::size_t nnodes = T_MESH_ELEMENT::n_all_nodes;

    //jacobian: J=[dx/dr dy/dr // dx/ds dy/ds]
    for (std::size_t k=0; k<nnodes; k++) {
        double const* const xyz = ele.getNode(k)->getCoords();
        // outer product of dNdr and xyz for a particular node
        for (std::size_t j_x=0; j_x<dim; j_x++) {
            for (std::size_t i_r=0; i_r<dim; i_r++) {
                shapemat.J(i_r,j_x) += shapemat.dNdr(i_r,k) * xyz[j_x];
            }
        }
    }

    shapemat.detJ = shapemat.J.determinant();
#ifndef NDEBUG
    if (shapemat.detJ<=.0)
        ERR("***error: det|J|=%e is not positive.\n", shapemat.detJ);
#endif
};

template <class T_MESH_ELEMENT, class T_SHAPE_FUNC, class T_SHAPE_MATRICES>
inline void computeMappingMatrices(
        const T_MESH_ELEMENT &ele,
        const double* natural_pt,
        T_SHAPE_MATRICES &shapemat,
        FieldType<ShapeMatrixType::N_J>)
{
    computeMappingMatrices<T_MESH_ELEMENT, T_SHAPE_FUNC, T_SHAPE_MATRICES>
        (ele, natural_pt, shapemat, FieldType<ShapeMatrixType::N>());
    computeMappingMatrices<T_MESH_ELEMENT, T_SHAPE_FUNC, T_SHAPE_MATRICES>
        (ele, natural_pt, shapemat, FieldType<ShapeMatrixType::DNDR_J>());
};

template <class T_MESH_ELEMENT, class T_SHAPE_FUNC, class T_SHAPE_MATRICES>
inline void computeMappingMatrices(
        const T_MESH_ELEMENT &ele,
        const double* natural_pt,
        T_SHAPE_MATRICES &shapemat,
        FieldType<ShapeMatrixType::DNDX>)
{
    computeMappingMatrices<T_MESH_ELEMENT, T_SHAPE_FUNC, T_SHAPE_MATRICES>
        (ele, natural_pt, shapemat, FieldType<ShapeMatrixType::DNDR_J>());

    if (shapemat.detJ>.0) {
        //J^-1, dshape/dx
        shapemat.invJ = shapemat.J.inverse();
        shapemat.dNdx = shapemat.invJ * shapemat.dNdr;
    }
};

template <class T_MESH_ELEMENT, class T_SHAPE_FUNC, class T_SHAPE_MATRICES>
inline void computeMappingMatrices(
        const T_MESH_ELEMENT &ele,
        const double* natural_pt,
        T_SHAPE_MATRICES &shapemat,
        FieldType<ShapeMatrixType::ALL>)
{
    computeMappingMatrices<T_MESH_ELEMENT, T_SHAPE_FUNC, T_SHAPE_MATRICES>
        (ele, natural_pt, shapemat, FieldType<ShapeMatrixType::N>());
    computeMappingMatrices<T_MESH_ELEMENT, T_SHAPE_FUNC, T_SHAPE_MATRICES>
        (ele, natural_pt, shapemat, FieldType<ShapeMatrixType::DNDX>());
};

} // detail

template <class T_MESH_ELEMENT, class T_SHAPE_FUNC, class T_SHAPE_MATRICES>
inline void NaturalCoordinatesMapping<
    T_MESH_ELEMENT,
    T_SHAPE_FUNC,
    T_SHAPE_MATRICES>
::computeShapeMatrices(
        const T_MESH_ELEMENT &ele,
        const double* natural_pt,
        T_SHAPE_MATRICES &shapemat)
{
    detail::computeMappingMatrices<
        T_MESH_ELEMENT,
        T_SHAPE_FUNC,
        T_SHAPE_MATRICES>
            (ele,
             natural_pt,
             shapemat,
             detail::FieldType<ShapeMatrixType::ALL>());
}

template <class T_MESH_ELEMENT, class T_SHAPE_FUNC, class T_SHAPE_MATRICES>
template <ShapeMatrixType T_SHAPE_MATRIX_TYPE>
inline void NaturalCoordinatesMapping<
    T_MESH_ELEMENT,
    T_SHAPE_FUNC,
    T_SHAPE_MATRICES>
::computeShapeMatrices(
        const T_MESH_ELEMENT &ele,
        const double* natural_pt,
        T_SHAPE_MATRICES &shapemat)
{
    detail::computeMappingMatrices<
        T_MESH_ELEMENT,
        T_SHAPE_FUNC,
        T_SHAPE_MATRICES>
            (ele,
             natural_pt,
             shapemat,
             detail::FieldType<T_SHAPE_MATRIX_TYPE>());
}

} // NumLib
