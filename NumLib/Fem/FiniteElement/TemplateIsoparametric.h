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


#ifndef TEMPLATEISOPARAMETRIC_H_
#define TEMPLATEISOPARAMETRIC_H_


#include <cassert>

#include "../CoordinatesMapping/ShapeMatrices.h"
#include "../CoordinatesMapping/NaturalCoordinatesMapping.h"

namespace NumLib
{

/**
 * \brief Template class for isoparametric elements
 *
 * \tparam ShapeFunctionType_   The shape function type.
 * \tparam ShapeMatrixTypes_    An aggregate of shape matrix types.
 */
template <
    class ShapeFunctionType_,
    class ShapeMatrixTypes_
    >
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
    using NaturalCoordsMappingType = NaturalCoordinatesMapping<
            MeshElementType, ShapeFunctionType, ShapeMatrices>;

    /**
     * Constructor without specifying a mesh element. setMeshElement() must be called afterwards.
     */
    TemplateIsoparametric()
    : _ele(nullptr)
    {
    }

    /**
     * Construct this object for the given mesh element.
     *
     * @param e                      Mesh element object
     */
    TemplateIsoparametric(const MeshElementType &e)
    : _ele(&e)
    {
    }

    ~TemplateIsoparametric() {}

    /// return current mesh element
    const MeshElementType* getMeshElement() const {return _ele;}

    /// Sets the mesh element
    void setMeshElement(const MeshElementType &e)
    {
        this->_ele = &e;
    }
    /**
     * compute shape functions
     *
     * @param natural_pt    position in natural coordinates
     * @param shape         evaluated shape function matrices
     * @param global_dim    global dimension
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

private:
    void computeIntegralMeasure(bool is_axially_symmetric,
                                ShapeMatrices& shape) const
    {
        if (!is_axially_symmetric) {
            shape.integralMeasure = 1.0;
            return;
        }
    }

    const MeshElementType* _ele;
};

} // NumLib

#endif //TEMPLATEISOPARAMETRIC_H_
