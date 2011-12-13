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

VtkCompositeGeoObjectFilter::VtkCompositeGeoObjectFilter( vtkAlgorithm* inputAlgorithm )
	: VtkCompositeFilter(inputAlgorithm), _type(GEOLIB::POINT), _threshold(vtkThreshold::New())
{
	if (inputAlgorithm->GetNumberOfInputPorts() && inputAlgorithm->GetNumberOfInputConnections(0))
	{
	  vtkAlgorithmOutput* ao = inputAlgorithm->GetInputConnection(0,0);

	  if (ao)
	  {
		vtkAlgorithm* parentAlg = ao->GetProducer();

	  if (dynamic_cast<VtkPolylinesSource*>(parentAlg) != NULL) 
		_type = GEOLIB::POLYLINE;
	  else if (dynamic_cast<VtkSurfacesSource*>(parentAlg) != NULL) 
		_type = GEOLIB::SURFACE;
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
	if (_type == GEOLIB::POINT)
	{
		 composite = new VtkCompositePointToGlyphFilter(surface);
		_outputAlgorithm = composite->GetOutputAlgorithm();
	}
	else if (_type == GEOLIB::POLYLINE)
	{
		composite = new VtkCompositeLineToTubeFilter(surface);
		_outputAlgorithm = composite->GetOutputAlgorithm();
	}
	else
		_outputAlgorithm = surface;

}

void VtkCompositeGeoObjectFilter::SetIndex(size_t idx) 
{
	_threshold->ThresholdBetween(idx, idx);
}