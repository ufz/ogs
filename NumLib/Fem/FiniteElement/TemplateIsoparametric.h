/**
 * \author Norihiro Watanabe
 * \date   2013-08-13
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
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
};

} // NumLib

#endif //TEMPLATEISOPARAMETRIC_H_
