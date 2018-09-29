/**
 * \author Norihiro Watanabe
 * \date   2013-08-13
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <boost/math/constants/constants.hpp>
#include <cassert>

#include "NumLib/Fem/CoordinatesMapping/NaturalCoordinatesMapping.h"
#include "NumLib/Fem/CoordinatesMapping/ShapeMatrices.h"

namespace NumLib
{
/**
 * \brief Template class for isoparametric elements
 *
 * \tparam ShapeFunctionType_   The shape function type.
 * \tparam ShapeMatrixTypes_    An aggregate of shape matrix types.
 */
template <class ShapeFunctionType_, class ShapeMatrixTypes_>
class TemplateIsoparametric
{
public:
    using ShapeFunctionType = ShapeFunctionType_;

    /// Coordinate mapping matrices type.
    using ShapeMatrices = typename ShapeMatrixTypes_::ShapeMatrices;

    /// Type of the underlying geometrical element.
    using MeshElementType = typename ShapeFunctionType_::MeshElement;

    /// Natural coordinates mapping tools specialization for specific
    /// MeshElement and ShapeFunction types.
    using NaturalCoordsMappingType =
        NaturalCoordinatesMapping<MeshElementType, ShapeFunctionType,
                                  ShapeMatrices>;

    /**
     * Constructor without specifying a mesh element. setMeshElement() must be
     * called afterwards.
     */
    TemplateIsoparametric() : _ele(nullptr) {}

    /**
     * Construct this object for the given mesh element.
     *
     * @param e                      Mesh element object
     */
    TemplateIsoparametric(const MeshElementType& e) : _ele(&e) {}
    ~TemplateIsoparametric() = default;

    /// return current mesh element
    const MeshElementType* getMeshElement() const { return _ele; }

    /// Sets the mesh element
    void setMeshElement(const MeshElementType& e) { this->_ele = &e; }
    /**
     * compute shape functions
     *
     * @param natural_pt            position in natural coordinates
     * @param shape                 evaluated shape function matrices
     * @param global_dim            global dimension
     * @param is_axially_symmetric  if true, the integral measure for cylinder
     *                              coordinates is used, and 1 otherwise.
     */
    void computeShapeFunctions(const double* natural_pt, ShapeMatrices& shape,
                               const unsigned global_dim,
                               bool is_axially_symmetric) const
    {
        NaturalCoordsMappingType::computeShapeMatrices(*_ele, natural_pt, shape,
                                                       global_dim);
        computeIntegralMeasure(is_axially_symmetric, shape);
    }

    /**
     * compute shape functions
     *
     * @tparam T_SHAPE_MATRIX_TYPE  shape matrix types to be calculated
     * @param natural_pt            position in natural coordinates
     * @param shape                 evaluated shape function matrices
     * @param global_dim            global dimension
     * @param is_axially_symmetric  if true, the integral measure for cylinder
     *                              coordinates is used, and 1 otherwise.
     */
    template <ShapeMatrixType T_SHAPE_MATRIX_TYPE>
    void computeShapeFunctions(const double* natural_pt, ShapeMatrices& shape,
                               const unsigned global_dim,
                               bool is_axially_symmetric) const
    {
        NaturalCoordsMappingType::template computeShapeMatrices<
            T_SHAPE_MATRIX_TYPE>(*_ele, natural_pt, shape, global_dim);
        computeIntegralMeasure(is_axially_symmetric, shape);
    }

    /// Interpolates the x coordinate of the element with the given shape
    /// matrix.
    double interpolateZerothCoordinate(
        typename ShapeMatrices::ShapeType const& N) const
    {
        auto* nodes = _ele->getNodes();
        typename ShapeMatrices::ShapeType rs(N.size());
        for (int i = 0; i < rs.size(); ++i)
        {
            rs[i] = (*nodes[i])[0];
        }
        auto const r = N.dot(rs);
        return r;
    }

    /// Interpolates the coordinates of the element with the given shape matrix.
    std::array<double, 3> interpolateCoordinates(
        typename ShapeMatrices::ShapeType const& N) const
    {
        auto* nodes = _ele->getNodes();
        typename ShapeMatrices::ShapeType rs(N.size());
        std::array<double, 3> interpolated_coords;
        for (int d = 0; d < 3; ++d)
        {
            for (int i = 0; i < rs.size(); ++i)
            {
                rs[i] = (*nodes[i])[d];
            }
            interpolated_coords[d] = N.dot(rs);
        }
        return interpolated_coords;
    }

private:
    void computeIntegralMeasure(bool is_axially_symmetric,
                                ShapeMatrices& shape) const
    {
        if (!is_axially_symmetric)
        {
            shape.integralMeasure = 1.0;
            return;
        }

        // Note: If an integration point is located at the rotation axis, r will
        // be zero which might lead to problems with the assembled equation
        // system.
        // E.g., for triangle elements, if an integration point is
        // located at edge of the unit triangle, it is possible that
        // r becomes zero.
        auto const r = interpolateZerothCoordinate(shape.N);
        shape.integralMeasure = boost::math::constants::two_pi<double>() * r;
    }

    const MeshElementType* _ele;
};

}  // namespace NumLib
