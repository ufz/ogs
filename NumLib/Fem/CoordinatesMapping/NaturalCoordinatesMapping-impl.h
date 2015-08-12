/**
 * \author Norihiro Watanabe
 * \date   2013-08-13
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */


#include <cassert>

#include <logog/include/logog.hpp>

#include "MeshLib/ElementCoordinatesMappingLocal.h"
#include "MeshLib/CoordinateSystem.h"

namespace NumLib
{

namespace detail
{

template <ShapeMatrixType FIELD_TYPE> struct FieldType {};

template <class T_MESH_ELEMENT, class T_SHAPE_FUNC, class T_SHAPE_MATRICES>
inline void computeMappingMatrices(
        const T_MESH_ELEMENT &/*ele*/,
        const double* natural_pt,
        const MeshLib::ElementCoordinatesMappingLocal &/*ele_local_coord*/,
        T_SHAPE_MATRICES &shapemat,
        FieldType<ShapeMatrixType::N>)
{
    T_SHAPE_FUNC::computeShapeFunction(natural_pt, shapemat.N);
}

template <class T_MESH_ELEMENT, class T_SHAPE_FUNC, class T_SHAPE_MATRICES>
inline void computeMappingMatrices(
        const T_MESH_ELEMENT &/*ele*/,
        const double* natural_pt,
        const MeshLib::ElementCoordinatesMappingLocal &/*ele_local_coord*/,
        T_SHAPE_MATRICES &shapemat,
        FieldType<ShapeMatrixType::DNDR>)
{
    double* const dNdr = shapemat.dNdr.data();
    T_SHAPE_FUNC::computeGradShapeFunction(natural_pt, dNdr);
}

template <class T_MESH_ELEMENT, class T_SHAPE_FUNC, class T_SHAPE_MATRICES>
inline void computeMappingMatrices(
        const T_MESH_ELEMENT &ele,
        const double* natural_pt,
        const MeshLib::ElementCoordinatesMappingLocal &ele_local_coord,
        T_SHAPE_MATRICES &shapemat,
        FieldType<ShapeMatrixType::DNDR_J>)
{
    computeMappingMatrices<T_MESH_ELEMENT, T_SHAPE_FUNC, T_SHAPE_MATRICES>
        (ele, natural_pt, ele_local_coord, shapemat, FieldType<ShapeMatrixType::DNDR>());

    const std::size_t dim = T_MESH_ELEMENT::dimension;
    const std::size_t nnodes = T_MESH_ELEMENT::n_all_nodes;

    //jacobian: J=[dx/dr dy/dr // dx/ds dy/ds]
    for (std::size_t k=0; k<nnodes; k++) {
        const MathLib::Point3d& mapped_pt = ele_local_coord.getMappedCoordinates(k);
        // outer product of dNdr and mapped_pt for a particular node
        for (std::size_t i_r=0; i_r<dim; i_r++) {
            for (std::size_t j_x=0; j_x<dim; j_x++) {
                shapemat.J(i_r,j_x) += shapemat.dNdr(i_r,k) * mapped_pt[j_x];
            }
        }
    }

    //shapemat.detJ = shapemat.J.determinant();
    if (ele.getDimension()==1)
        shapemat.detJ = shapemat.J(0,0);
    else if (ele.getDimension()==2)
        shapemat.detJ = shapemat.J(0,0) * shapemat.J(1,1) - shapemat.J(0,1) * shapemat.J(1,0);
    else if (ele.getDimension()==3)
        shapemat.detJ = shapemat.J(0,0) * (shapemat.J(1,1) * shapemat.J(2,2) - shapemat.J(2,1) * shapemat.J(1,2))
                 + shapemat.J(2,0) * (shapemat.J(0,1) * shapemat.J(1,2) - shapemat.J(1,1) * shapemat.J(0,2))
                 + shapemat.J(1,0) * (shapemat.J(0,2) * shapemat.J(2,1) - shapemat.J(2,2) * shapemat.J(0,1));

#ifndef NDEBUG
    if (shapemat.detJ<=.0)
        ERR("***error: det|J|=%e is not positive.\n", shapemat.detJ);
#endif
}

template <class T_MESH_ELEMENT, class T_SHAPE_FUNC, class T_SHAPE_MATRICES>
inline void computeMappingMatrices(
        const T_MESH_ELEMENT &ele,
        const double* natural_pt,
        const MeshLib::ElementCoordinatesMappingLocal &ele_local_coord,
        T_SHAPE_MATRICES &shapemat,
        FieldType<ShapeMatrixType::N_J>)
{
    computeMappingMatrices<T_MESH_ELEMENT, T_SHAPE_FUNC, T_SHAPE_MATRICES>
        (ele, natural_pt, ele_local_coord, shapemat, FieldType<ShapeMatrixType::N>());
    computeMappingMatrices<T_MESH_ELEMENT, T_SHAPE_FUNC, T_SHAPE_MATRICES>
        (ele, natural_pt, ele_local_coord, shapemat, FieldType<ShapeMatrixType::DNDR_J>());
}

template <class T_MESH_ELEMENT, class T_SHAPE_FUNC, class T_SHAPE_MATRICES>
inline void computeMappingMatrices(
        const T_MESH_ELEMENT &ele,
        const double* natural_pt,
        const MeshLib::ElementCoordinatesMappingLocal &ele_local_coord,
        T_SHAPE_MATRICES &shapemat,
        FieldType<ShapeMatrixType::DNDX>)
{
    computeMappingMatrices<T_MESH_ELEMENT, T_SHAPE_FUNC, T_SHAPE_MATRICES>
        (ele, natural_pt, ele_local_coord, shapemat, FieldType<ShapeMatrixType::DNDR_J>());

    if (shapemat.detJ>.0) {
        //J^-1, dshape/dx
        //shapemat.invJ.noalias() = shapemat.J.inverse();
        if (ele.getDimension()==1) {
            shapemat.invJ(0,0) = 1./shapemat.detJ;
        } else if (ele.getDimension()==2) {
            shapemat.invJ(0,0) = shapemat.J(1,1);
            shapemat.invJ(0,1) = -shapemat.J(0,1);
            shapemat.invJ(1,0) = -shapemat.J(1,0);
            shapemat.invJ(1,1) = shapemat.J(0,0);
            shapemat.invJ *= 1./shapemat.detJ;
        } else if (ele.getDimension()==3) {
            shapemat.invJ(0,0) =  shapemat.J(1,1) * shapemat.J(2,2) - shapemat.J(2,1) * shapemat.J(1,2);
            shapemat.invJ(0,1) =  shapemat.J(0,2) * shapemat.J(2,1) - shapemat.J(0,1) * shapemat.J(2,2);
            shapemat.invJ(0,2) =  shapemat.J(0,1) * shapemat.J(1,2) - shapemat.J(0,2) * shapemat.J(1,1);
            //
            shapemat.invJ(1,0) =  shapemat.J(1,2) * shapemat.J(2,0) - shapemat.J(2,2) * shapemat.J(1,0);
            shapemat.invJ(1,1) =  shapemat.J(0,0) * shapemat.J(2,2) - shapemat.J(2,0) * shapemat.J(0,2);
            shapemat.invJ(1,2) =  shapemat.J(0,2) * shapemat.J(1,0) - shapemat.J(1,2) * shapemat.J(0,0);
            //
            shapemat.invJ(2,0) =  shapemat.J(1,0) * shapemat.J(2,1) - shapemat.J(2,0) * shapemat.J(1,1);
            shapemat.invJ(2,1) =  shapemat.J(0,1) * shapemat.J(2,0) - shapemat.J(2,1) * shapemat.J(0,0);
            shapemat.invJ(2,2) =  shapemat.J(0,0) * shapemat.J(1,1) - shapemat.J(1,0) * shapemat.J(0,1);
            shapemat.invJ /= shapemat.detJ;
        }

        auto const nnodes(shapemat.dNdr.cols());
        auto const ele_dim(shapemat.dNdr.rows());
        assert(shapemat.dNdr.rows()==ele.getDimension());
        const unsigned global_dim(ele_local_coord.getGlobalCoordinateSystem().getDimension());
        if (global_dim==ele_dim) {
            shapemat.dNdx.topLeftCorner(ele_dim, nnodes).noalias() = shapemat.invJ * shapemat.dNdr;
        } else {
            auto const& matR = ele_local_coord.getRotationMatrixToGlobal(); // 3 x 3
            auto invJ_dNdr = shapemat.invJ * shapemat.dNdr;
            auto dshape_global = matR.topLeftCorner(3u, ele_dim) * invJ_dNdr; //3 x nnodes
            shapemat.dNdx = dshape_global.topLeftCorner(global_dim, nnodes);;
        }
    }
}

template <class T_MESH_ELEMENT, class T_SHAPE_FUNC, class T_SHAPE_MATRICES>
inline void computeMappingMatrices(
        const T_MESH_ELEMENT &ele,
        const double* natural_pt,
        const MeshLib::ElementCoordinatesMappingLocal &ele_local_coord,
        T_SHAPE_MATRICES &shapemat,
        FieldType<ShapeMatrixType::ALL>)
{
    computeMappingMatrices<T_MESH_ELEMENT, T_SHAPE_FUNC, T_SHAPE_MATRICES>
        (ele, natural_pt, ele_local_coord, shapemat, FieldType<ShapeMatrixType::N>());
    computeMappingMatrices<T_MESH_ELEMENT, T_SHAPE_FUNC, T_SHAPE_MATRICES>
        (ele, natural_pt, ele_local_coord, shapemat, FieldType<ShapeMatrixType::DNDX>());
}

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
    const MeshLib::CoordinateSystem coords(ele);
    const MeshLib::ElementCoordinatesMappingLocal ele_local_coord(ele, coords);

    detail::computeMappingMatrices<
        T_MESH_ELEMENT,
        T_SHAPE_FUNC,
        T_SHAPE_MATRICES>
            (ele,
             natural_pt,
             ele_local_coord,
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
    const MeshLib::CoordinateSystem coords(ele);
    const MeshLib::ElementCoordinatesMappingLocal ele_local_coord(ele, coords);

    detail::computeMappingMatrices<
        T_MESH_ELEMENT,
        T_SHAPE_FUNC,
        T_SHAPE_MATRICES>
            (ele,
             natural_pt,
             ele_local_coord,
             shapemat,
             detail::FieldType<T_SHAPE_MATRIX_TYPE>());
}

} // NumLib
