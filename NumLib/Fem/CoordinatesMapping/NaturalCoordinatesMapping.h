/**
 * \author Norihiro Watanabe
 * \date   2013-08-13
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */


#ifndef NATURALCOORDINATESMAPPING_H_
#define NATURALCOORDINATESMAPPING_H_


#include "ShapeMatrices.h"

namespace NumLib
{
namespace detail
{
template <class T_MESH_ELEMENT,
          class T_SHAPE_FUNC,
          class T_SHAPE_MATRICES,
          ShapeMatrixType T_SHAPE_MATRIX_TYPE>
void naturalCoordinatesMappingComputeShapeMatrices(const T_MESH_ELEMENT& ele,
                                                   const double* natural_pt,
                                                   T_SHAPE_MATRICES& shapemat);
}  // namespace detail

/**
 * Coordinates mapping tools for natural coordinates
 *
 * This class also supports coordinates mapping of mixed dimensional elements,
 * e.g. line elements in 2D space. Details of the mapping method can be found in
 * \cite Kolditz2001 .
 *
 * @tparam T_MESH_ELEMENT       Mesh element type
 * @tparam T_SHAPE_FUNC         Shape function class
 * @tparam T_SHAPE_MATRICES     Shape matrices class
 */
template <class T_MESH_ELEMENT, class T_SHAPE_FUNC, class T_SHAPE_MATRICES>
struct NaturalCoordinatesMapping
{
    /**
     * compute all mapping matrices at the given location in natural coordinates
     *
     * @param ele               Mesh element object
     * @param natural_pt        Location in natural coordinates (r,s,t)
     * @param shapemat          Shape matrix data where calculated shape functions are stored
     */
    static void computeShapeMatrices(const T_MESH_ELEMENT& ele,
                                     const double* natural_pt,
                                     T_SHAPE_MATRICES& shapemat)
    {
        computeShapeMatrices<ShapeMatrixType::ALL>(ele, natural_pt, shapemat);
    }

    /**
     * compute specified mapping matrices at the given location in natural
     * coordinates
     *
     * @tparam T_SHAPE_MATRIX_TYPE  Mapping matrix types to be calculated
     * @param ele                   Mesh element object
     * @param natural_pt            Location in natural coordinates (r,s,t)
     * @param shapemat              Shape matrix data where calculated shape
     * functions are stored
     */
    template <ShapeMatrixType T_SHAPE_MATRIX_TYPE>
    static void computeShapeMatrices(const T_MESH_ELEMENT& ele,
                                     const double* natural_pt,
                                     T_SHAPE_MATRICES& shapemat)
    {
        detail::naturalCoordinatesMappingComputeShapeMatrices<
            T_MESH_ELEMENT,
            T_SHAPE_FUNC,
            T_SHAPE_MATRICES,
            T_SHAPE_MATRIX_TYPE>(ele, natural_pt, shapemat);
    }
};

} // NumLib

// #include "NaturalCoordinatesMapping-impl.h"


#endif //NATURALCOORDINATESMAPPING_H_
