/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "NaturalCoordinatesMapping.h"

#include <cassert>

#include "BaseLib/Error.h"

#include "MeshLib/ElementCoordinatesMappingLocal.h"
#include "MeshLib/Elements/TemplateElement.h"
#include "MeshLib/Elements/HexRule20.h"
#include "MeshLib/Elements/HexRule8.h"
#include "MeshLib/Elements/LineRule2.h"
#include "MeshLib/Elements/LineRule3.h"
#include "MeshLib/Elements/PointRule1.h"
#include "MeshLib/Elements/PrismRule15.h"
#include "MeshLib/Elements/PrismRule6.h"
#include "MeshLib/Elements/PyramidRule13.h"
#include "MeshLib/Elements/PyramidRule5.h"
#include "MeshLib/Elements/QuadRule4.h"
#include "MeshLib/Elements/QuadRule8.h"
#include "MeshLib/Elements/QuadRule9.h"
#include "MeshLib/Elements/TetRule10.h"
#include "MeshLib/Elements/TetRule4.h"
#include "MeshLib/Elements/TriRule3.h"
#include "MeshLib/Elements/TriRule6.h"

#include "NumLib/Fem/ShapeFunction/ShapePoint1.h"
#include "NumLib/Fem/ShapeFunction/ShapeLine2.h"
#include "NumLib/Fem/ShapeFunction/ShapeLine3.h"
#include "NumLib/Fem/ShapeFunction/ShapeTri3.h"
#include "NumLib/Fem/ShapeFunction/ShapeTri6.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad4.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad8.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad9.h"
#include "NumLib/Fem/ShapeFunction/ShapeHex8.h"
#include "NumLib/Fem/ShapeFunction/ShapeHex20.h"
#include "NumLib/Fem/ShapeFunction/ShapeTet4.h"
#include "NumLib/Fem/ShapeFunction/ShapeTet10.h"
#include "NumLib/Fem/ShapeFunction/ShapePrism6.h"
#include "NumLib/Fem/ShapeFunction/ShapePrism15.h"
#include "NumLib/Fem/ShapeFunction/ShapePyra5.h"
#include "NumLib/Fem/ShapeFunction/ShapePyra13.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"

#include "ShapeMatrices.h"

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
inline
typename std::enable_if<T_SHAPE_FUNC::DIM!=0>::type
computeMappingMatrices(
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
inline
typename std::enable_if<T_SHAPE_FUNC::DIM==0>::type
computeMappingMatrices(
        const T_MESH_ELEMENT &/*ele*/,
        const double* /*natural_pt*/,
        const MeshLib::ElementCoordinatesMappingLocal &/*ele_local_coord*/,
        T_SHAPE_MATRICES &/*shapemat*/,
        FieldType<ShapeMatrixType::DNDR>)
{
}

template <class T_MESH_ELEMENT, class T_SHAPE_FUNC, class T_SHAPE_MATRICES>
inline
typename std::enable_if<T_SHAPE_FUNC::DIM!=0>::type
computeMappingMatrices(
        const T_MESH_ELEMENT &ele,
        const double* natural_pt,
        const MeshLib::ElementCoordinatesMappingLocal &ele_local_coord,
        T_SHAPE_MATRICES &shapemat,
        FieldType<ShapeMatrixType::DNDR_J>)
{
    computeMappingMatrices<T_MESH_ELEMENT, T_SHAPE_FUNC, T_SHAPE_MATRICES>
        (ele, natural_pt, ele_local_coord, shapemat, FieldType<ShapeMatrixType::DNDR>());

    auto const dim = T_MESH_ELEMENT::dimension;
    auto const nnodes = T_MESH_ELEMENT::n_all_nodes;

    //jacobian: J=[dx/dr dy/dr // dx/ds dy/ds]
    for (auto k = decltype(nnodes){0}; k<nnodes; k++) {
        const MathLib::Point3d& mapped_pt = ele_local_coord.getMappedCoordinates(k);
        // outer product of dNdr and mapped_pt for a particular node
        for (auto i_r = decltype(dim){0}; i_r<dim; i_r++) {
            for (auto j_x = decltype(dim){0}; j_x<dim; j_x++) {
                shapemat.J(i_r,j_x) += shapemat.dNdr(i_r,k) * mapped_pt[j_x];
            }
        }
    }

    shapemat.detJ = shapemat.J.determinant();

    if (shapemat.detJ<=.0)
        OGS_FATAL("det J = %e is not positive.\n", shapemat.detJ);
}
template <class T_MESH_ELEMENT, class T_SHAPE_FUNC, class T_SHAPE_MATRICES>
inline
typename std::enable_if<T_SHAPE_FUNC::DIM==0>::type
computeMappingMatrices(
        const T_MESH_ELEMENT &/*ele*/,
        const double* /*natural_pt*/,
        const MeshLib::ElementCoordinatesMappingLocal &/*ele_local_coord*/,
        T_SHAPE_MATRICES &shapemat,
        FieldType<ShapeMatrixType::DNDR_J>)
{
    shapemat.detJ = 1.0;
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
inline
typename std::enable_if<T_SHAPE_FUNC::DIM!=0>::type
computeMappingMatrices(
        const T_MESH_ELEMENT &ele,
        const double* natural_pt,
        const MeshLib::ElementCoordinatesMappingLocal &ele_local_coord,
        T_SHAPE_MATRICES &shapemat,
        FieldType<ShapeMatrixType::DNDX>)
{
    computeMappingMatrices<T_MESH_ELEMENT, T_SHAPE_FUNC, T_SHAPE_MATRICES>
        (ele, natural_pt, ele_local_coord, shapemat, FieldType<ShapeMatrixType::DNDR_J>());

    if (shapemat.detJ > 0) {
        //J^-1, dshape/dx
        shapemat.invJ.noalias() = shapemat.J.inverse();

        auto const nnodes(shapemat.dNdr.cols());
        auto const ele_dim(shapemat.dNdr.rows());
        assert(shapemat.dNdr.rows()==ele.getDimension());
        const unsigned global_dim = ele_local_coord.getGlobalDimension();
        if (global_dim==ele_dim) {
            shapemat.dNdx.topLeftCorner(ele_dim, nnodes).noalias() = shapemat.invJ * shapemat.dNdr;
        } else {
            auto const& matR = ele_local_coord.getRotationMatrixToGlobal(); // 3 x 3
            auto invJ_dNdr = shapemat.invJ * shapemat.dNdr;
            auto dshape_global = matR.topLeftCorner(3u, ele_dim) * invJ_dNdr; //3 x nnodes
            shapemat.dNdx = dshape_global.topLeftCorner(global_dim, nnodes);;
        }
    } else {
        OGS_FATAL("det J = %e is not positive.\n", shapemat.detJ);
    }
}

template <class T_MESH_ELEMENT, class T_SHAPE_FUNC, class T_SHAPE_MATRICES>
inline
typename std::enable_if<T_SHAPE_FUNC::DIM==0>::type
computeMappingMatrices(
       const T_MESH_ELEMENT &ele,
       const double* natural_pt,
       const MeshLib::ElementCoordinatesMappingLocal &ele_local_coord,
       T_SHAPE_MATRICES &shapemat,
       FieldType<ShapeMatrixType::DNDX>)
{
    computeMappingMatrices<T_MESH_ELEMENT, T_SHAPE_FUNC, T_SHAPE_MATRICES>
        (ele, natural_pt, ele_local_coord, shapemat, FieldType<ShapeMatrixType::DNDR_J>());
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

template <class T_MESH_ELEMENT,
          class T_SHAPE_FUNC,
          class T_SHAPE_MATRICES,
          ShapeMatrixType T_SHAPE_MATRIX_TYPE>
void naturalCoordinatesMappingComputeShapeMatrices(const T_MESH_ELEMENT& ele,
                                                   const double* natural_pt,
                                                   T_SHAPE_MATRICES& shapemat,
                                                   const unsigned global_dim)
{
    const MeshLib::ElementCoordinatesMappingLocal ele_local_coord(ele, global_dim);

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

#define OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_PART(        \
    RULE, SHAPE, DIM, WHICHPART, SHAPEMATRIXPOLICY)              \
    template void naturalCoordinatesMappingComputeShapeMatrices< \
        MeshLib::TemplateElement<MeshLib::RULE>,                 \
        NumLib::SHAPE,                                           \
        SHAPEMATRIXPOLICY<NumLib::SHAPE, DIM>::ShapeMatrices,    \
        ShapeMatrixType::WHICHPART>(                             \
        MeshLib::TemplateElement<MeshLib::RULE> const&,          \
        double const*,                                           \
        SHAPEMATRIXPOLICY<NumLib::SHAPE, DIM>::ShapeMatrices&,   \
        const unsigned global_dim)

#define OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(RULE, SHAPE) \
    OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_PART(                \
        RULE, SHAPE, 0, ALL, EigenDynamicShapeMatrixPolicy);         \
    /* Those instantiations are needed in unit tests only */         \
    OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_PART(                \
        RULE, SHAPE, 0, N, EigenDynamicShapeMatrixPolicy);           \
    OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_PART(                \
        RULE, SHAPE, 0, DNDR, EigenDynamicShapeMatrixPolicy);        \
    OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_PART(                \
        RULE, SHAPE, 0, N_J, EigenDynamicShapeMatrixPolicy);         \
    OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_PART(                \
        RULE, SHAPE, 0, DNDR_J, EigenDynamicShapeMatrixPolicy);      \
    OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_PART(                \
        RULE, SHAPE, 0, DNDX, EigenDynamicShapeMatrixPolicy)

#define OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(RULE, SHAPE, DIM) \
    OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_PART(                     \
        RULE, SHAPE, DIM, ALL, EigenFixedShapeMatrixPolicy);              \
    /* Those instantiations are needed in unit tests only */              \
    OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_PART(                     \
        RULE, SHAPE, DIM, N, EigenFixedShapeMatrixPolicy);                \
    OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_PART(                     \
        RULE, SHAPE, DIM, DNDR, EigenFixedShapeMatrixPolicy);             \
    OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_PART(                     \
        RULE, SHAPE, DIM, N_J, EigenFixedShapeMatrixPolicy);              \
    OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_PART(                     \
        RULE, SHAPE, DIM, DNDR_J, EigenFixedShapeMatrixPolicy);           \
    OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_PART(                     \
        RULE, SHAPE, DIM, DNDX, EigenFixedShapeMatrixPolicy)

OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(HexRule20, ShapeHex20);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(HexRule8, ShapeHex8);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(LineRule2, ShapeLine2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(LineRule3, ShapeLine3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(PointRule1, ShapePoint1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(PrismRule15, ShapePrism15);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(PrismRule6, ShapePrism6);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(PyramidRule13, ShapePyra13);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(PyramidRule5, ShapePyra5);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(QuadRule4, ShapeQuad4);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(QuadRule8, ShapeQuad8);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(QuadRule9, ShapeQuad9);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(TetRule10, ShapeTet10);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(TetRule4, ShapeTet4);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(TriRule3, ShapeTri3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(TriRule6, ShapeTri6);

OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(HexRule20, ShapeHex20, 1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(HexRule8, ShapeHex8, 1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(LineRule2, ShapeLine2, 1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(LineRule3, ShapeLine3, 1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(PointRule1, ShapePoint1, 1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(PrismRule15, ShapePrism15, 1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(PrismRule6, ShapePrism6, 1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(PyramidRule13, ShapePyra13, 1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(PyramidRule5, ShapePyra5, 1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(QuadRule4, ShapeQuad4, 1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(QuadRule8, ShapeQuad8, 1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(QuadRule9, ShapeQuad9, 1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(TetRule10, ShapeTet10, 1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(TetRule4, ShapeTet4, 1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(TriRule3, ShapeTri3, 1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(TriRule6, ShapeTri6, 1);

OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(HexRule20, ShapeHex20, 2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(HexRule8, ShapeHex8, 2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(LineRule2, ShapeLine2, 2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(LineRule3, ShapeLine3, 2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(PointRule1, ShapePoint1, 2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(PrismRule15, ShapePrism15, 2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(PrismRule6, ShapePrism6, 2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(PyramidRule13, ShapePyra13, 2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(PyramidRule5, ShapePyra5, 2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(QuadRule4, ShapeQuad4, 2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(QuadRule8, ShapeQuad8, 2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(QuadRule9, ShapeQuad9, 2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(TetRule10, ShapeTet10, 2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(TetRule4, ShapeTet4, 2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(TriRule3, ShapeTri3, 2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(TriRule6, ShapeTri6, 2);

OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(HexRule20, ShapeHex20, 3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(HexRule8, ShapeHex8, 3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(LineRule2, ShapeLine2, 3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(LineRule3, ShapeLine3, 3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(PointRule1, ShapePoint1, 3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(PrismRule15, ShapePrism15, 3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(PrismRule6, ShapePrism6, 3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(PyramidRule13, ShapePyra13, 3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(PyramidRule5, ShapePyra5, 3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(QuadRule4, ShapeQuad4, 3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(QuadRule8, ShapeQuad8, 3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(QuadRule9, ShapeQuad9, 3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(TetRule10, ShapeTet10, 3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(TetRule4, ShapeTet4, 3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(TriRule3, ShapeTri3, 3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(TriRule6, ShapeTri6, 3);

} // detail

} // NumLib
