/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "NaturalCoordinatesMapping.h"

#include <cassert>
#ifndef NDEBUG
#include <iostream>
#endif  // NDEBUG

#include "BaseLib/Error.h"
#include "MeshLib/ElementCoordinatesMappingLocal.h"
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
#include "MeshLib/Elements/TemplateElement.h"
#include "MeshLib/Elements/TetRule10.h"
#include "MeshLib/Elements/TetRule4.h"
#include "MeshLib/Elements/TriRule3.h"
#include "MeshLib/Elements/TriRule6.h"
#include "NumLib/Fem/ShapeFunction/ShapeHex20.h"
#include "NumLib/Fem/ShapeFunction/ShapeHex8.h"
#include "NumLib/Fem/ShapeFunction/ShapeLine2.h"
#include "NumLib/Fem/ShapeFunction/ShapeLine3.h"
#include "NumLib/Fem/ShapeFunction/ShapePoint1.h"
#include "NumLib/Fem/ShapeFunction/ShapePrism15.h"
#include "NumLib/Fem/ShapeFunction/ShapePrism6.h"
#include "NumLib/Fem/ShapeFunction/ShapePyra13.h"
#include "NumLib/Fem/ShapeFunction/ShapePyra5.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad4.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad8.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad9.h"
#include "NumLib/Fem/ShapeFunction/ShapeTet10.h"
#include "NumLib/Fem/ShapeFunction/ShapeTet4.h"
#include "NumLib/Fem/ShapeFunction/ShapeTri3.h"
#include "NumLib/Fem/ShapeFunction/ShapeTri6.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ShapeMatrices.h"

namespace NumLib
{
namespace detail
{
template <ShapeMatrixType FIELD_TYPE>
struct FieldType
{
};

template <class T_SHAPE_FUNC, class T_SHAPE_MATRICES>
inline void computeMappingMatrices(
    const MeshLib::Element& /*ele*/,
    const double* natural_pt,
    const MeshLib::ElementCoordinatesMappingLocal& /*ele_local_coord*/,
    T_SHAPE_MATRICES& shapemat,
    FieldType<ShapeMatrixType::N> /*unused*/)
{
    T_SHAPE_FUNC::computeShapeFunction(natural_pt, shapemat.N);
}

template <class T_SHAPE_FUNC, class T_SHAPE_MATRICES>
inline void computeMappingMatrices(
    const MeshLib::Element& /*ele*/,
    [[maybe_unused]] const double* natural_pt,
    const MeshLib::ElementCoordinatesMappingLocal& /*ele_local_coord*/,
    [[maybe_unused]] T_SHAPE_MATRICES& shapemat,
    FieldType<ShapeMatrixType::DNDR> /*unused*/)
{
    if constexpr (T_SHAPE_FUNC::DIM != 0)
    {
        double* const dNdr = shapemat.dNdr.data();
        T_SHAPE_FUNC::computeGradShapeFunction(natural_pt, dNdr);
    }
}

static void checkJacobianDeterminant(const double detJ,
                                     MeshLib::Element const& element)
{
    if (detJ > 0)
    {  // The usual case
        return;
    }

    if (detJ < 0)
    {
        ERR("det J = {:g} is negative for element {:d}.",
            detJ,
            element.getID());
#ifndef NDEBUG
        std::cerr << element << "\n";
#endif  // NDEBUG
        OGS_FATAL(
            "Please check whether the node numbering of the element is correct,"
            "or additional elements (like boundary elements) are still present "
            "in the mesh.");
    }

    if (detJ == 0)
    {
        ERR("det J is zero for element {:d}.", element.getID());
#ifndef NDEBUG
        std::cerr << element << "\n";
#endif  // NDEBUG
        OGS_FATAL(
            "Please check whether:\n"
            "\t the element nodes may have the same coordinates,\n"
            "\t or the coordinates of all nodes are not given on the x-axis "
            "for a 1D problem or in the xy-plane in the 2D case.\n"
            "The first case can occur, if not all boundary elements"
            "were removed from the bulk mesh.");
    }
}

template <class T_SHAPE_FUNC, class T_SHAPE_MATRICES>
inline void computeMappingMatrices(
    const MeshLib::Element& ele,
    [[maybe_unused]] const double* natural_pt,
    const MeshLib::ElementCoordinatesMappingLocal& ele_local_coord,
    T_SHAPE_MATRICES& shapemat,
    FieldType<ShapeMatrixType::DNDR_J> /*unused*/)
{
    if constexpr (T_SHAPE_FUNC::DIM != 0)
    {
        computeMappingMatrices<T_SHAPE_FUNC, T_SHAPE_MATRICES>(
            ele,
            natural_pt,
            ele_local_coord,
            shapemat,
            FieldType<ShapeMatrixType::DNDR>());

        auto const dim = T_SHAPE_FUNC::DIM;
        auto const nnodes = T_SHAPE_FUNC::NPOINTS;

        // jacobian: J=[dx/dr dy/dr // dx/ds dy/ds]
        for (auto k = decltype(nnodes){0}; k < nnodes; k++)
        {
            const MathLib::Point3d& mapped_pt =
                ele_local_coord.getMappedCoordinates(k);
            // outer product of dNdr and mapped_pt for a particular node
            for (auto i_r = decltype(dim){0}; i_r < dim; i_r++)
            {
                for (auto j_x = decltype(dim){0}; j_x < dim; j_x++)
                {
                    shapemat.J(i_r, j_x) +=
                        shapemat.dNdr(i_r, k) * mapped_pt[j_x];
                }
            }
        }

        shapemat.detJ = shapemat.J.determinant();
        checkJacobianDeterminant(shapemat.detJ, ele);
    }
    else
    {
        shapemat.detJ = 1.0;
    }
}

template <class T_SHAPE_FUNC, class T_SHAPE_MATRICES>
inline void computeMappingMatrices(
    const MeshLib::Element& ele,
    const double* natural_pt,
    const MeshLib::ElementCoordinatesMappingLocal& ele_local_coord,
    T_SHAPE_MATRICES& shapemat,
    FieldType<ShapeMatrixType::N_J> /*unused*/)
{
    computeMappingMatrices<T_SHAPE_FUNC, T_SHAPE_MATRICES>(
        ele,
        natural_pt,
        ele_local_coord,
        shapemat,
        FieldType<ShapeMatrixType::N>());
    computeMappingMatrices<T_SHAPE_FUNC, T_SHAPE_MATRICES>(
        ele,
        natural_pt,
        ele_local_coord,
        shapemat,
        FieldType<ShapeMatrixType::DNDR_J>());
}

template <class T_SHAPE_FUNC, class T_SHAPE_MATRICES>
inline void computeMappingMatrices(
    const MeshLib::Element& ele,
    const double* natural_pt,
    const MeshLib::ElementCoordinatesMappingLocal& ele_local_coord,
    T_SHAPE_MATRICES& shapemat,
    FieldType<ShapeMatrixType::DNDX> /*unused*/)
{
    computeMappingMatrices<T_SHAPE_FUNC, T_SHAPE_MATRICES>(
        ele,
        natural_pt,
        ele_local_coord,
        shapemat,
        FieldType<ShapeMatrixType::DNDR_J>());
    if constexpr (T_SHAPE_FUNC::DIM != 0)
    {
        checkJacobianDeterminant(shapemat.detJ, ele);

        // J^-1, dshape/dx
        shapemat.invJ.noalias() = shapemat.J.inverse();

        assert(shapemat.dNdr.rows() == ele.getDimension());
        const unsigned global_dim = ele_local_coord.getGlobalDimension();
        if (global_dim == T_SHAPE_FUNC::DIM)
        {
            shapemat.dNdx
                .template topLeftCorner<T_SHAPE_FUNC::DIM,
                                        T_SHAPE_FUNC::NPOINTS>()
                .noalias() = shapemat.invJ * shapemat.dNdr;
        }
        else
        {
            auto const& matR =
                (ele_local_coord.getRotationMatrixToGlobal().topLeftCorner(
                     global_dim, T_SHAPE_FUNC::DIM))
                    .eval();

            auto const invJ_dNdr = shapemat.invJ * shapemat.dNdr;
            auto const dshape_global = matR * invJ_dNdr;
            shapemat.dNdx =
                dshape_global.topLeftCorner(global_dim, T_SHAPE_FUNC::NPOINTS);
        }
    }
}

template <class T_SHAPE_FUNC, class T_SHAPE_MATRICES>
inline void computeMappingMatrices(
    const MeshLib::Element& ele,
    const double* natural_pt,
    const MeshLib::ElementCoordinatesMappingLocal& ele_local_coord,
    T_SHAPE_MATRICES& shapemat,
    FieldType<ShapeMatrixType::ALL> /*unused*/)
{
    computeMappingMatrices<T_SHAPE_FUNC, T_SHAPE_MATRICES>(
        ele,
        natural_pt,
        ele_local_coord,
        shapemat,
        FieldType<ShapeMatrixType::N>());
    computeMappingMatrices<T_SHAPE_FUNC, T_SHAPE_MATRICES>(
        ele,
        natural_pt,
        ele_local_coord,
        shapemat,
        FieldType<ShapeMatrixType::DNDX>());
}

template <class T_SHAPE_FUNC,
          class T_SHAPE_MATRICES,
          ShapeMatrixType T_SHAPE_MATRIX_TYPE>
void naturalCoordinatesMappingComputeShapeMatrices(const MeshLib::Element& ele,
                                                   const double* natural_pt,
                                                   T_SHAPE_MATRICES& shapemat,
                                                   const unsigned global_dim)
{
    const MeshLib::ElementCoordinatesMappingLocal ele_local_coord(ele,
                                                                  global_dim);

    detail::computeMappingMatrices<T_SHAPE_FUNC, T_SHAPE_MATRICES>(
        ele,
        natural_pt,
        ele_local_coord,
        shapemat,
        detail::FieldType<T_SHAPE_MATRIX_TYPE>());
}

#define OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_PART(        \
    SHAPE, DIM, WHICHPART, SHAPEMATRIXPOLICY)                    \
    template void naturalCoordinatesMappingComputeShapeMatrices< \
        NumLib::SHAPE,                                           \
        SHAPEMATRIXPOLICY<NumLib::SHAPE, DIM>::ShapeMatrices,    \
        ShapeMatrixType::WHICHPART>(                             \
        MeshLib::Element const&,                                 \
        double const*,                                           \
        SHAPEMATRIXPOLICY<NumLib::SHAPE, DIM>::ShapeMatrices&,   \
        const unsigned global_dim)

#define OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(SHAPE) \
    OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_PART(          \
        SHAPE, 0, ALL, EigenDynamicShapeMatrixPolicy);         \
    /* Those instantiations are needed in unit tests only */   \
    OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_PART(          \
        SHAPE, 0, N, EigenDynamicShapeMatrixPolicy);           \
    OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_PART(          \
        SHAPE, 0, DNDR, EigenDynamicShapeMatrixPolicy);        \
    OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_PART(          \
        SHAPE, 0, N_J, EigenDynamicShapeMatrixPolicy);         \
    OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_PART(          \
        SHAPE, 0, DNDR_J, EigenDynamicShapeMatrixPolicy);      \
    OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_PART(          \
        SHAPE, 0, DNDX, EigenDynamicShapeMatrixPolicy)

#define OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(SHAPE, DIM) \
    OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_PART(               \
        SHAPE, DIM, ALL, EigenFixedShapeMatrixPolicy);              \
    /* Those instantiations are needed in unit tests only */        \
    OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_PART(               \
        SHAPE, DIM, N, EigenFixedShapeMatrixPolicy);                \
    OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_PART(               \
        SHAPE, DIM, DNDR, EigenFixedShapeMatrixPolicy);             \
    OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_PART(               \
        SHAPE, DIM, N_J, EigenFixedShapeMatrixPolicy);              \
    OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_PART(               \
        SHAPE, DIM, DNDR_J, EigenFixedShapeMatrixPolicy);           \
    OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_PART(               \
        SHAPE, DIM, DNDX, EigenFixedShapeMatrixPolicy)

OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(ShapeHex20);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(ShapeHex8);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(ShapeLine2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(ShapeLine3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(ShapePoint1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(ShapePrism15);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(ShapePrism6);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(ShapePyra13);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(ShapePyra5);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(ShapeQuad4);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(ShapeQuad8);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(ShapeQuad9);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(ShapeTet10);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(ShapeTet4);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(ShapeTri3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_DYN(ShapeTri6);

OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeHex20, 1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeHex8, 1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeLine2, 1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeLine3, 1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapePoint1, 1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapePrism15, 1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapePrism6, 1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapePyra13, 1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapePyra5, 1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeQuad4, 1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeQuad8, 1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeQuad9, 1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeTet10, 1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeTet4, 1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeTri3, 1);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeTri6, 1);

OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeHex20, 2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeHex8, 2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeLine2, 2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeLine3, 2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapePoint1, 2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapePrism15, 2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapePrism6, 2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapePyra13, 2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapePyra5, 2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeQuad4, 2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeQuad8, 2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeQuad9, 2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeTet10, 2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeTet4, 2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeTri3, 2);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeTri6, 2);

OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeHex20, 3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeHex8, 3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeLine2, 3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeLine3, 3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapePoint1, 3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapePrism15, 3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapePrism6, 3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapePyra13, 3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapePyra5, 3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeQuad4, 3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeQuad8, 3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeQuad9, 3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeTet10, 3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeTet4, 3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeTri3, 3);
OGS_INSTANTIATE_NATURAL_COORDINATES_MAPPING_FIX(ShapeTri6, 3);

}  // namespace detail

}  // namespace NumLib
