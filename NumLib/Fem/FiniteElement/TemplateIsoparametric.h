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
 * \tparam T_ELEMENT            Mesh element class
 * \tparam T_SHAPE              Shape function
 * \tparam T_INTEGRAL           Integration method
 * \tparam T_NODAL_VECTOR       Nodal vector class
 * \tparam T_DIM_NODAL_MATRIX   Matrix class for a size of dim * nnodes
 * \tparam T_DIM_MATRIX         Matrix class for a size of dim * dim
 */
template <
    class T_MESH_ELEMENT,
    class T_SHAPE,
    class T_INTEGRAL,
    class T_NODAL_VECTOR,
    class T_DIM_NODAL_MATRIX,
    class T_DIM_MATRIX
    >
class TemplateIsoparametric
{
public:
    typedef T_MESH_ELEMENT MeshElementType;
    typedef T_SHAPE ShapeFunctionType;
    typedef T_INTEGRAL IntegrationMethod;
    typedef T_NODAL_VECTOR NodalVectorType;
    typedef T_DIM_NODAL_MATRIX DimNodalMatrixType;
    typedef T_DIM_MATRIX DimMatrixType;
    typedef ShapeMatrices<NodalVectorType, DimNodalMatrixType, DimMatrixType> ShapeMatricesType;
    typedef NaturalCoordinatesMapping<MeshElementType, ShapeFunctionType, ShapeMatricesType> NaturalCoordsMappingType;

    /**
     * Constructor without specifying a mesh element. resetMeshElement() must be called afterwards.
     *
     * @param integration_order      default 2.
     */
    explicit TemplateIsoparametric(std::size_t integration_order=2)
    : _ele(nullptr)
    {
        this->_integration.setIntegrationOrder(integration_order);
    }

    /**
     * Construct this object for the given mesh element.
     *
     * @param e                      Mesh element object
     * @param integration_order      default 2.
     */
    TemplateIsoparametric(const MeshElementType &e, std::size_t integration_order=2)
    : _ele(&e)
    {
        this->_integration.setIntegrationOrder(integration_order);
    }

    ///
    ~TemplateIsoparametric() {}

    /// return current mesh element
    const MeshElementType* getMeshElement() const {return _ele;}

    /// reset a mesh element
    void resetMeshElement(const MeshElementType &e)
    {
        this->_ele = &e;
    }

    /// return an integration method
    const IntegrationMethod& getIntegrationMethod() const {return _integration;}

    /**
     * compute shape functions
     *
     * @param natural_pt    position in natural coordinates
     * @param shape         evaluated shape function matrices
     */
    void computeShapeFunctions(const double *natural_pt, ShapeMatricesType &shape) const
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
    void computeShapeFunctions(const double *natural_pt, ShapeMatricesType &shape) const
    {
        NaturalCoordsMappingType::template computeShapeMatrices<T_SHAPE_MATRIX_TYPE>(*_ele, natural_pt, shape);
    }


private:
    const MeshElementType* _ele;
    IntegrationMethod _integration;
};

} // NumLib

#endif //TEMPLATEISOPARAMETRIC_H_
