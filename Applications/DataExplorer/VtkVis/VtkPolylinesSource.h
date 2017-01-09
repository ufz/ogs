/**
 * \file
 * \author Lars Bilke
 * \date   2010-02-02
 * \brief  Definition of the VtkPolylinesSource class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef VTKPOLYLINESSOURCE_H
#define VTKPOLYLINESSOURCE_H

// ** INCLUDES **
#include "VtkAlgorithmProperties.h"
#include <vtkPolyDataAlgorithm.h>

// forward declaration
namespace GeoLib {
class Polyline;
}


/**
 * \brief VtkPolylinesSource is a VTK source object for the visualisation of
 * polyline data. As a vtkPolyDataAlgorithm it outputs polygonal data.
 */
class VtkPolylinesSource : public vtkPolyDataAlgorithm, public VtkAlgorithmProperties
{
public:
    /// Create new objects with New() because of VTKs object reference counting.
    static VtkPolylinesSource* New();

    vtkTypeMacro(VtkPolylinesSource,vtkPolyDataAlgorithm);

    /// Sets the polyline vector.
    void setPolylines(const std::vector<GeoLib::Polyline*>* polylines) { _polylines = polylines; }

    /// Prints its data on a stream.
    void PrintSelf(ostream& os, vtkIndent indent) override;

    virtual void SetUserProperty(QString name, QVariant value);

protected:
    VtkPolylinesSource();
    ~VtkPolylinesSource();

    /// Computes the polygonal data object.
    int RequestData(vtkInformation* request,
                    vtkInformationVector** inputVector,
                    vtkInformationVector* outputVector) override;

    int RequestInformation(vtkInformation* request,
                           vtkInformationVector** inputVector,
                           vtkInformationVector* outputVector) override;

    /// The polylines to visualize.
    const std::vector<GeoLib::Polyline*>* _polylines;

private:
};

#endif // VTKPOLYLINESSOURCE_H
