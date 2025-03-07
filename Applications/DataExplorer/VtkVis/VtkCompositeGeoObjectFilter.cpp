/**
 * \file
 * \author Karsten Rink
 * \date   2011-12-02
 * \brief  Implementation of the VtkCompositeGeoObjectFilter class.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkCompositeGeoObjectFilter.h"

#include <vtkAlgorithmOutput.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkThreshold.h>

#include "VtkCompositeLineToTubeFilter.h"
#include "VtkCompositePointToGlyphFilter.h"
#include "VtkPolylinesSource.h"
#include "VtkStationSource.h"
#include "VtkSurfacesSource.h"

VtkCompositeGeoObjectFilter::VtkCompositeGeoObjectFilter(
    vtkAlgorithm* inputAlgorithm)
    : VtkCompositeFilter(inputAlgorithm),
      _type(GeoLib::GEOTYPE::POINT),
      _threshold(vtkThreshold::New())
{
    if (inputAlgorithm->GetNumberOfInputPorts() &&
        inputAlgorithm->GetNumberOfInputConnections(0))
    {
        vtkAlgorithmOutput* ao = inputAlgorithm->GetInputConnection(0, 0);

        if (ao)
        {
            vtkAlgorithm* parentAlg = ao->GetProducer();

            if (dynamic_cast<VtkPolylinesSource*>(parentAlg) != nullptr)
            {
                _type = GeoLib::GEOTYPE::POLYLINE;
            }
            else if (dynamic_cast<VtkSurfacesSource*>(parentAlg) != nullptr)
            {
                _type = GeoLib::GEOTYPE::SURFACE;
            }
            else if (dynamic_cast<VtkStationSource*>(parentAlg) != nullptr)
            {
                /* TODO
                if (dynamic_cast<VtkStationSource*>(parentAlg)->getType() ==
                GeoLib::Station::StationType::BOREHOLE) _type =
                GeoLib::GEOTYPE::POLYLINE;
                */
            }
        }
        this->init();
    }
}

VtkCompositeGeoObjectFilter::~VtkCompositeGeoObjectFilter() = default;

void VtkCompositeGeoObjectFilter::init()
{
    this->_inputDataObjectType = VTK_POLY_DATA;
    this->_outputDataObjectType = VTK_POLY_DATA;

    _threshold->SetInputConnection(_inputAlgorithm->GetOutputPort());
    _threshold->SetSelectedComponent(0);
    _threshold->SetThresholdFunction(
        vtkThreshold::ThresholdType::THRESHOLD_BETWEEN);
    _threshold->SetLowerThreshold(0);
    _threshold->SetUpperThreshold(0);

    vtkDataSetSurfaceFilter* surface = vtkDataSetSurfaceFilter::New();
    surface->SetInputConnection(_threshold->GetOutputPort());

    VtkCompositeFilter* composite;
    if (_type == GeoLib::GEOTYPE::POINT)
    {
        composite = new VtkCompositePointToGlyphFilter(surface);
        composite->SetUserProperty("Radius", this->GetInitialRadius());
        _outputAlgorithm = composite->GetOutputAlgorithm();
    }
    else if (_type == GeoLib::GEOTYPE::POLYLINE)
    {
        composite = new VtkCompositeLineToTubeFilter(surface);
        composite->SetUserProperty("Radius", this->GetInitialRadius());
        _outputAlgorithm = composite->GetOutputAlgorithm();
    }
    else
    {
        _outputAlgorithm = surface;
    }
}

void VtkCompositeGeoObjectFilter::SetIndex(std::size_t idx)
{
    double const d_idx = static_cast<double>(idx);
    _threshold->SetThresholdFunction(
        vtkThreshold::ThresholdType::THRESHOLD_BETWEEN);
    _threshold->SetLowerThreshold(d_idx);
    _threshold->SetUpperThreshold(d_idx);
}
