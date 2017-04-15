/**
 * \file
 * \author Lars Bilke
 * \date   2010-02-03
 * \brief  Implementation of the VtkPointsSource class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "VtkPointsSource.h"

#include <logog/include/logog.hpp>

#include <vtkCellArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkCellData.h>
#include <vtkProperty.h>

#include "Applications/DataHolderLib/Color.h"

vtkStandardNewMacro(VtkPointsSource);

VtkPointsSource::VtkPointsSource() : _points(nullptr)
{
    _removable = false; // From VtkAlgorithmProperties
    this->SetNumberOfInputPorts(0);

    const DataHolderLib::Color c = DataHolderLib::getRandomColor();
    GetProperties()->SetColor(c[0] / 255.0,c[1] / 255.0,c[2] / 255.0);
}

void VtkPointsSource::PrintSelf( ostream& os, vtkIndent indent )
{
    this->Superclass::PrintSelf(os,indent);

    if (_points->size() == 0)
        return;

    os << indent << "== VtkPointsSource ==" << "\n";

    int i = 0;
    for (std::vector<GeoLib::Point*>::const_iterator it = _points->begin();
         it != _points->end(); ++it)
    {
        const double* coords = (*it)->getCoords();
        os << indent << "Point " << i << " (" << coords[0] << ", " << coords[1] << ", " <<
        coords[2] << ")\n";
        i++;
    }
}

int VtkPointsSource::RequestData( vtkInformation* request,
                                  vtkInformationVector** inputVector,
                                  vtkInformationVector* outputVector )
{
    (void)request;
    (void)inputVector;

    if (!_points)
        return 0;
    int numPoints = _points->size();
    if (numPoints == 0)
    {
        ERR("VtkPointsSource::RequestData(): Size of point vector is 0");
        return 0;
    }

    vtkSmartPointer<vtkInformation> outInfo = outputVector->GetInformationObject(0);
    vtkSmartPointer<vtkPolyData> output =
            vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> newVerts = vtkSmartPointer<vtkCellArray>::New();
    newPoints->SetNumberOfPoints(numPoints);
    newVerts->Allocate(numPoints);

    vtkSmartPointer<vtkIntArray> pointIDs = vtkSmartPointer<vtkIntArray>::New();
    pointIDs->SetNumberOfComponents(1);
    pointIDs->SetNumberOfValues(numPoints);
    pointIDs->SetName("PointIDs");

    if (outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()) > 0)
        return 1;

    // Generate points and vertices
    unsigned i = 0;
    for (std::vector<GeoLib::Point*>::const_iterator it = _points->begin();
         it != _points->end(); ++it)
    {
        double coords[3] = {(*(*it))[0], (*(*it))[1], (*(*it))[2]};
        newPoints->SetPoint(i, coords);
        newVerts->InsertNextCell(1);
        newVerts->InsertCellPoint(i);

        pointIDs->SetValue(i, i);
        i++;
    }

    output->SetPoints(newPoints);
    output->SetVerts(newVerts);
    output->GetCellData()->AddArray(pointIDs);
    output->GetCellData()->SetActiveAttribute("PointIDs", vtkDataSetAttributes::SCALARS);

    return 1;
}

int VtkPointsSource::RequestInformation( vtkInformation* /*request*/,
                                         vtkInformationVector** /*inputVector*/,
                                         vtkInformationVector* /*outputVector*/ )
{
    return 1;
}

void VtkPointsSource::SetUserProperty( QString name, QVariant value )
{
    Q_UNUSED(name);
    Q_UNUSED(value);
}
