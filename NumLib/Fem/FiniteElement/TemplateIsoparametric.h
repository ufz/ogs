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


#ifndef TEMPLATEISOPARAMETRIC_H_
#define TEMPLATEISOPARAMETRIC_H_


#include <vector>
#include <cassert>

#include "NumLib/Fem/CoordinatesMapping/NaturalCoordinatesMapping.h"

namespace NumLib
{

/**
 * \brief Template class for isoparametric elements
 *
 * \tparam T_ELEMENT            Mesh element class
 * \tparam T_SHAPE              Shape function
 * \tparam T_INTEGRAL           Integration method
 * \tparam T_EXTRAPOLATE        Extrapolation method
 * \tparam T_NODAL_VECTOR       Nodal vector class
 * \tparam T_DIM_NODAL_MATRIX   Matrix class for a size of dim * nnodes
 * \tparam T_DIM_MATRIX         Matrix class for a size of dim * dim
 */
template <
    class T_MESH_ELEMENT,
    class T_SHAPE,
    class T_INTEGRAL,
    class T_EXTRAPOLATE,
    class T_NODAL_VECTOR,
    class T_DIM_NODAL_MATRIX,
    class T_DIM_MATRIX
    >
class TemplateIsoparametric
{
public:
    typedef T_MESH_ELEMENT MeshElementType;
    typedef T_SHAPE ShapeFunctionType;
    typedef T_INTEGRAL IntegrationType;
    typedef T_EXTRAPOLATE ExtrapolationType;
    typedef T_NODAL_VECTOR NodalVectorType;
    typedef T_DIM_NODAL_MATRIX DimNodalMatrixType;
    typedef T_DIM_MATRIX DimMatrixType;
    typedef ShapeData<NodalVectorType, DimNodalMatrixType, DimMatrixType> ShapeDataType;

    /**
     * Constructor without specifying a mesh element. resetMeshElement() must be called afterwards.
     *
     * @param integration_points_level      the sampling level (default 2)
     */
    explicit TemplateIsoparametric(std::size_t integration_points_level=2)
    : _ele(nullptr)
    {
        this->_integration.setSamplingLevel(integration_points_level);
    }

    /**
     * Construct this object for the given mesh element.
     *
     * @param e                             Mesh element object
     * @param integration_points_level      the sampling level (default 2)
     */
    TemplateIsoparametric(const MeshElementType &e, std::size_t integration_points_level=2)
    : _ele(&e)
    {
        this->resetMeshElement(e);
        this->_integration.setSamplingLevel(integration_points_level);
    }

    ///
    ~TemplateIsoparametric() {}

    /// return current mesh element
    const MeshElementType* getMeshElement() const {return _ele;}

    /// reset a mesh element
    void resetMeshElement(const MeshElementType &e)
    {
        this->_ele = &e;
        this->_mapping.reset(e);
    }

    /// return an integration method
    const IntegrationType& getIntegrationMethod() const {return _integration;}

    /**
     * compute shape functions
     *
     * @param x         point in natural coordinates
     * @param shape     evaluated shape data
     */
    void computeShapeFunctions(const double *x, ShapeDataType &shape) const
    {
        _mapping.computeMappingMatrices(x, shape);
    }

    /// make interpolation from nodal values
    double interpolate(const ShapeDataType& shape, NodalVectorType &nodal_values) const
    {
        assert(nodal_values.size()==this->_ele->getNNodes());
        return shape.N.transpose().dot(nodal_values);
    }

    /// extrapolate integration point values to nodal values
    template <class T_VEC>
    void extrapolate(const T_VEC &gp_values, NodalVectorType &nodal_values) const
    {
        ExtrapolationType::extrapolate(gp_values, nodal_values);
    }

private:
    const MeshElementType* _ele;
    NaturalCoordinatesMapping<MeshElementType, ShapeFunctionType, ShapeDataType> _mapping;
    IntegrationType _integration;
};

} //end

#endif //TEMPLATEISOPARAMETRIC_H_
