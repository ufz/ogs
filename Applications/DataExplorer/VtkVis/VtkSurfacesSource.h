/**
 * \file
 * \author Lars Bilke
 * \date   2010-02-03
 * \brief  Definition of the VtkSurfacesSource class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

// ** INCLUDES **
#include "VtkAlgorithmProperties.h"
#include <vtkPolyDataAlgorithm.h>

#include "Surface.h"

/**
 * \brief VTK source object for the visualisation of surfaces.
 * Technically, surfaces are displayed as triangulated polydata objects.
 */
class VtkSurfacesSource : public vtkPolyDataAlgorithm, public VtkAlgorithmProperties
{
public:
    /// Create new objects with New() because of VTKs object reference counting.
    static VtkSurfacesSource* New();

    vtkTypeMacro(VtkSurfacesSource,vtkPolyDataAlgorithm);

    /// Sets the surfaces vector
    void setSurfaces(const std::vector<GeoLib::Surface*>* surfaces) { _surfaces = surfaces; }

    /// Prints its data on a stream.
    void PrintSelf(ostream& os, vtkIndent indent) override;

    /**
     * \brief Generates random colors for each surface.
     */
    //ogsUserPropertyMacro(ColorBySurface,bool);

    void SetUserProperty(QString name, QVariant value) override;

protected:
    VtkSurfacesSource();
    ~VtkSurfacesSource() override = default;

    /// Computes the polygonal data object.
    int RequestData(vtkInformation* request,
                    vtkInformationVector** inputVector,
                    vtkInformationVector* outputVector) override;

    int RequestInformation(vtkInformation* request,
                           vtkInformationVector** inputVector,
                           vtkInformationVector* outputVector) override;

    /// The surfaces to visualize
    const std::vector<GeoLib::Surface*>* _surfaces;

private:
};
