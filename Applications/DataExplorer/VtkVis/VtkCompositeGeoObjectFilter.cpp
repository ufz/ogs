/**
 * \file
 * \author Karsten Rink
 * \date   2011-12-02
 * \brief  Implementation of the VtkCompositeGeoObjectFilter class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkCompositeGeoObjectFilter.h"

#include <vtkDataSetSurfaceFilter.h>
#include <vtkSmartPointer.h>
#include <vtkThreshold.h>
#include <vtkAlgorithmOutput.h>

#include "VtkPolylinesSource.h"
#include "VtkSurfacesSource.h"
#include "VtkStationSource.h"
#include "VtkCompositePointToGlyphFilter.h"
#include "VtkCompositeLineToTubeFilter.h"

#include <vtkPointData.h>

VtkCompositeGeoObjectFilter::VtkCompositeGeoObjectFilter( vtkAlgorithm* inputAlgorithm )
    : VtkCompositeFilter(inputAlgorithm), type_(GeoLib::GEOTYPE::POINT), threshold_(vtkThreshold::New())
{
    if (inputAlgorithm->GetNumberOfInputPorts() && inputAlgorithm->GetNumberOfInputConnections(0))
    {
      vtkAlgorithmOutput* ao = inputAlgorithm->GetInputConnection(0,0);

      if (ao)
      {
        vtkAlgorithm* parentAlg = ao->GetProducer();

        if (dynamic_cast<VtkPolylinesSource*>(parentAlg) != nullptr)
        {
            type_ = GeoLib::GEOTYPE::POLYLINE;
        }
        else if (dynamic_cast<VtkSurfacesSource*>(parentAlg) != nullptr)
        {
            type_ = GeoLib::GEOTYPE::SURFACE;
        }
        else if (dynamic_cast<VtkStationSource*>(parentAlg) != nullptr)
        {
            /* TODO
            if (dynamic_cast<VtkStationSource*>(parentAlg)->getType() == GeoLib::Station::StationType::BOREHOLE)
                type_ = GeoLib::GEOTYPE::POLYLINE;
            */
        }
      }
      this->init();
    }
}

VtkCompositeGeoObjectFilter::~VtkCompositeGeoObjectFilter() = default;

void VtkCompositeGeoObjectFilter::init()
{
    this->inputDataObjectType_ = VTK_POLY_DATA;
    this->outputDataObjectType_ = VTK_POLY_DATA;

    threshold_->SetInputConnection(inputAlgorithm_->GetOutputPort());
    threshold_->SetSelectedComponent(0);
    threshold_->ThresholdBetween(0,0);

    vtkDataSetSurfaceFilter* surface = vtkDataSetSurfaceFilter::New();
    surface->SetInputConnection(threshold_->GetOutputPort());

    VtkCompositeFilter* composite;
    if (type_ == GeoLib::GEOTYPE::POINT)
    {
        composite = new VtkCompositePointToGlyphFilter(surface);
        composite->SetUserProperty("Radius", this->GetInitialRadius());
        outputAlgorithm_ = composite->GetOutputAlgorithm();
    }
    else if (type_ == GeoLib::GEOTYPE::POLYLINE)
    {
        composite = new VtkCompositeLineToTubeFilter(surface);
        composite->SetUserProperty("Radius", this->GetInitialRadius());
        outputAlgorithm_ = composite->GetOutputAlgorithm();
    }
    else
    {
        outputAlgorithm_ = surface;
    }
}

void VtkCompositeGeoObjectFilter::SetIndex(std::size_t idx)
{
    double const d_idx = static_cast<double>(idx);
    threshold_->ThresholdBetween(d_idx, d_idx);
}


