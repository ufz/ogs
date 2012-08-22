/**
 * \file VtkCompositeGeoObjectFilter.cpp
 * 2011/12/02 KR Initial implementation
 *
 * Implementation of VtkCompositeGeoObjectFilter class
 */

// ** INCLUDES **
#include "VtkCompositeGeoObjectFilter.h"

#include <vtkDataSetSurfaceFilter.h>
#include <vtkSmartPointer.h>
#include <vtkThreshold.h>
#include <vtkAlgorithmOutput.h>

#include "VtkPolylinesSource.h"
#include "VtkSurfacesSource.h"
#include "VtkCompositePointToGlyphFilter.h"
#include "VtkCompositeLineToTubeFilter.h"

#include <vtkPointData.h>

VtkCompositeGeoObjectFilter::VtkCompositeGeoObjectFilter( vtkAlgorithm* inputAlgorithm )
	: VtkCompositeFilter(inputAlgorithm), _type(GeoLib::POINT), _threshold(vtkThreshold::New())
{
	if (inputAlgorithm->GetNumberOfInputPorts() && inputAlgorithm->GetNumberOfInputConnections(0))
	{
	  vtkAlgorithmOutput* ao = inputAlgorithm->GetInputConnection(0,0);

	  if (ao)
	  {
		vtkAlgorithm* parentAlg = ao->GetProducer();

	  if (dynamic_cast<VtkPolylinesSource*>(parentAlg) != NULL)
		_type = GeoLib::POLYLINE;
	  else if (dynamic_cast<VtkSurfacesSource*>(parentAlg) != NULL)
		_type = GeoLib::SURFACE;
	  }

	}

	this->init();
}

VtkCompositeGeoObjectFilter::~VtkCompositeGeoObjectFilter()
{
}

void VtkCompositeGeoObjectFilter::init()
{
	this->_inputDataObjectType = VTK_POLY_DATA;
	this->_outputDataObjectType = VTK_POLY_DATA;

	_threshold->SetInputConnection(_inputAlgorithm->GetOutputPort());
	_threshold->SetSelectedComponent(0);
	_threshold->ThresholdBetween(0,0);

	vtkDataSetSurfaceFilter* surface = vtkDataSetSurfaceFilter::New();
	surface->SetInputConnection(_threshold->GetOutputPort());

	VtkCompositeFilter* composite;
	if (_type == GeoLib::POINT)
	{
		composite = new VtkCompositePointToGlyphFilter(surface);
		composite->SetUserProperty("Radius", this->GetInitialRadius());
		_outputAlgorithm = composite->GetOutputAlgorithm();
	}
	else if (_type == GeoLib::POLYLINE)
	{
		composite = new VtkCompositeLineToTubeFilter(surface);
		composite->SetUserProperty("Radius", this->GetInitialRadius());
		_outputAlgorithm = composite->GetOutputAlgorithm();

	}
	else
		_outputAlgorithm = surface;

}

void VtkCompositeGeoObjectFilter::SetIndex(size_t idx)
{
	_threshold->ThresholdBetween(idx, idx);
}

float VtkCompositeGeoObjectFilter::GetInitialRadius() const
{
	double bounding_box[6];
	static_cast<vtkPolyData*>(this->_inputAlgorithm->GetOutputDataObject(0))->GetBounds(bounding_box);
	double x_diff = fabs(bounding_box[0]-bounding_box[1]);
	double y_diff = fabs(bounding_box[2]-bounding_box[3]);
	double z_diff = fabs(bounding_box[4]-bounding_box[5]);

	double max = (x_diff > y_diff) ? x_diff : y_diff;
	max = (max > z_diff) ? max : z_diff;

	return max/200.0;
}
