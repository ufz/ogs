/**
 * \file
 * \author Lars Bilke
 * \date   2010-10-19
 * \brief  Definition of the VtkCompositeFilter class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "VtkAlgorithmProperties.h"

class vtkAlgorithm;

/**
 * @brief Is used to combine several filter in one VtkVisPipelineItem. You can
 * use vtk filter and custom filter. To subclass this you have to implement the
 * init() function. There you combine the filters. Make sure to set the members
 * _inputDataObjectType, _outputDataObjectType and _outputAlgorithm. Make also
 * sure to implement VtkAlgorithmProperties::SetUserProperty() and
 * VtkAlgorithmProperties::SetUserVectorProperty().
 *
 * Allocate vtk objects inside init() with vtkSmartPointer except for the last
 * filter. This filter must also be set to _outputAlgorithm, e.g.
 *
 *     MyVtkFilter* lastFilter = MyVtkFilter::New();
 *     ...(do something here)
 *     _outputAlgorithm = lastFilter;
 *
 * Create user properties with `ogsUserPropertyMacro` or `ogsUserVecxPropertyMacro`
 * and initialize these properties inside the constructor with
 * `this->Set[Property Name](value)`.
 * See VtkCompositeThresholdFilter for an example.
 */
class VtkCompositeFilter : public VtkAlgorithmProperties
{
public:
    /// @brief Constructor.
    /// @param inputAlgorithm The algorithm to attach this filter to.
    VtkCompositeFilter(vtkAlgorithm* inputAlgorithm);

    /// @brief Destructor.
    ~VtkCompositeFilter() override;

    /// @return the type of the data input.
    /// Can be compared with
    /// - VTK_POLY_DATA
    /// - VTK_STRUCTURED_POINTS
    /// - VTK_STRUCTURED_GRID
    /// - VTK_RECTILINEAR_GRID
    /// - VTK_UNSTRUCTURED_GRID
    /// - VTK_IMAGE_DATA
    /// - VTK_DATA_SET
    int GetInputDataObjectType() const { return _inputDataObjectType; }

    /// @return the type of the data output.
    int GetOutputDataObjectType() const { return _outputDataObjectType; }

    /// @return the last algorithm in this composite filter.
    vtkAlgorithm* GetOutputAlgorithm() const { return _outputAlgorithm; }

protected:
    /// Calculates a 1/200th of the largest extension of the bounding box (this is used as default radius for various filters)
    float GetInitialRadius() const;

    /// See [vtkSetGet.h](https://github.com/Kitware/VTK/blob/master/Common/Core/vtkSetGet.h)
    /// for the defines
    int _inputDataObjectType;
    int _outputDataObjectType;

    vtkAlgorithm* _inputAlgorithm;
    vtkAlgorithm* _outputAlgorithm;

    virtual void init() = 0;
};
