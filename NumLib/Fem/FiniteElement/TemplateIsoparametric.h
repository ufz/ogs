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
     */
    void computeShapeFunctions(const double *natural_pt, ShapeMatrices &shape) const
    {
        NaturalCoordsMappingType::computeShapeMatrices(*_ele, natural_pt, shape);
    }

    /**
     * compute shape functions
     *
     * @tparam T_SHAPE_MATRIX_TYPE  shape matrix types to be calculated
     * @param natural_pt            position in natural coordinates
     * @param shape                 evaluated shape function matrices
     */
    template <ShapeMatrixType T_SHAPE_MATRIX_TYPE>
    void computeShapeFunctions(const double *natural_pt, ShapeMatrices &shape) const
    {
        NaturalCoordsMappingType::template computeShapeMatrices<T_SHAPE_MATRIX_TYPE>(*_ele, natural_pt, shape);
    }


private:
    const MeshElementType* _ele;
};

} // NumLib

#include "NumLib/Fem/ShapeMatrixPolicy.h"

extern template class NumLib::TemplateIsoparametric<NumLib::ShapeHex20  , EigenFixedShapeMatrixPolicy<NumLib::ShapeHex20  , 3>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeHex8   , EigenFixedShapeMatrixPolicy<NumLib::ShapeHex8   , 3>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeLine2  , EigenFixedShapeMatrixPolicy<NumLib::ShapeLine2  , 3>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeLine3  , EigenFixedShapeMatrixPolicy<NumLib::ShapeLine3  , 3>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapePrism15, EigenFixedShapeMatrixPolicy<NumLib::ShapePrism15, 3>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapePrism6 , EigenFixedShapeMatrixPolicy<NumLib::ShapePrism6 , 3>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapePyra13 , EigenFixedShapeMatrixPolicy<NumLib::ShapePyra13 , 3>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapePyra5  , EigenFixedShapeMatrixPolicy<NumLib::ShapePyra5  , 3>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeQuad4  , EigenFixedShapeMatrixPolicy<NumLib::ShapeQuad4  , 3>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeQuad8  , EigenFixedShapeMatrixPolicy<NumLib::ShapeQuad8  , 3>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeQuad9  , EigenFixedShapeMatrixPolicy<NumLib::ShapeQuad9  , 3>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeTet10  , EigenFixedShapeMatrixPolicy<NumLib::ShapeTet10  , 3>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeTet4   , EigenFixedShapeMatrixPolicy<NumLib::ShapeTet4   , 3>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeTri3   , EigenFixedShapeMatrixPolicy<NumLib::ShapeTri3   , 3>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeTri6   , EigenFixedShapeMatrixPolicy<NumLib::ShapeTri6   , 3>>;

extern template class NumLib::TemplateIsoparametric<NumLib::ShapeLine2  , EigenFixedShapeMatrixPolicy<NumLib::ShapeLine2  , 2>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeLine3  , EigenFixedShapeMatrixPolicy<NumLib::ShapeLine3  , 2>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeQuad4  , EigenFixedShapeMatrixPolicy<NumLib::ShapeQuad4  , 2>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeQuad8  , EigenFixedShapeMatrixPolicy<NumLib::ShapeQuad8  , 2>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeQuad9  , EigenFixedShapeMatrixPolicy<NumLib::ShapeQuad9  , 2>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeTri3   , EigenFixedShapeMatrixPolicy<NumLib::ShapeTri3   , 2>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeTri6   , EigenFixedShapeMatrixPolicy<NumLib::ShapeTri6   , 2>>;

extern template class NumLib::TemplateIsoparametric<NumLib::ShapeLine2  , EigenFixedShapeMatrixPolicy<NumLib::ShapeLine2  , 1>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeLine3  , EigenFixedShapeMatrixPolicy<NumLib::ShapeLine3  , 1>>;

extern template class NumLib::TemplateIsoparametric<NumLib::ShapeHex20  , EigenDynamicShapeMatrixPolicy<NumLib::ShapeHex20  , 3>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeHex8   , EigenDynamicShapeMatrixPolicy<NumLib::ShapeHex8   , 3>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeLine2  , EigenDynamicShapeMatrixPolicy<NumLib::ShapeLine2  , 3>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeLine3  , EigenDynamicShapeMatrixPolicy<NumLib::ShapeLine3  , 3>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapePrism15, EigenDynamicShapeMatrixPolicy<NumLib::ShapePrism15, 3>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapePrism6 , EigenDynamicShapeMatrixPolicy<NumLib::ShapePrism6 , 3>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapePyra13 , EigenDynamicShapeMatrixPolicy<NumLib::ShapePyra13 , 3>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapePyra5  , EigenDynamicShapeMatrixPolicy<NumLib::ShapePyra5  , 3>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeQuad4  , EigenDynamicShapeMatrixPolicy<NumLib::ShapeQuad4  , 3>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeQuad8  , EigenDynamicShapeMatrixPolicy<NumLib::ShapeQuad8  , 3>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeQuad9  , EigenDynamicShapeMatrixPolicy<NumLib::ShapeQuad9  , 3>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeTet10  , EigenDynamicShapeMatrixPolicy<NumLib::ShapeTet10  , 3>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeTet4   , EigenDynamicShapeMatrixPolicy<NumLib::ShapeTet4   , 3>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeTri3   , EigenDynamicShapeMatrixPolicy<NumLib::ShapeTri3   , 3>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeTri6   , EigenDynamicShapeMatrixPolicy<NumLib::ShapeTri6   , 3>>;

extern template class NumLib::TemplateIsoparametric<NumLib::ShapeLine2  , EigenDynamicShapeMatrixPolicy<NumLib::ShapeLine2  , 2>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeLine3  , EigenDynamicShapeMatrixPolicy<NumLib::ShapeLine3  , 2>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeQuad4  , EigenDynamicShapeMatrixPolicy<NumLib::ShapeQuad4  , 2>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeQuad8  , EigenDynamicShapeMatrixPolicy<NumLib::ShapeQuad8  , 2>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeQuad9  , EigenDynamicShapeMatrixPolicy<NumLib::ShapeQuad9  , 2>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeTri3   , EigenDynamicShapeMatrixPolicy<NumLib::ShapeTri3   , 2>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeTri6   , EigenDynamicShapeMatrixPolicy<NumLib::ShapeTri6   , 2>>;

extern template class NumLib::TemplateIsoparametric<NumLib::ShapeLine2  , EigenDynamicShapeMatrixPolicy<NumLib::ShapeLine2  , 1>>;
extern template class NumLib::TemplateIsoparametric<NumLib::ShapeLine3  , EigenDynamicShapeMatrixPolicy<NumLib::ShapeLine3  , 1>>;

#endif //TEMPLATEISOPARAMETRIC_H_
