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

#include "VtkCompositePointToGlyphFilter.h"
#include "VtkCompositeLineToTubeFilter.h"

VtkCompositeGeoObjectFilter::VtkCompositeGeoObjectFilter( vtkAlgorithm* inputAlgorithm )
	: VtkCompositeFilter(inputAlgorithm), _index(0), _type(GEOLIB::INVALID)
{
	this->init();
}

VtkCompositeGeoObjectFilter::~VtkCompositeGeoObjectFilter()
{
}

void VtkCompositeGeoObjectFilter::init()
{
	this->_inputDataObjectType = VTK_POLY_DATA;
	this->_outputDataObjectType = VTK_POLY_DATA;

	vtkThreshold* threshold = vtkThreshold::New();
	threshold->SetInputConnection(_inputAlgorithm->GetOutputPort());
	threshold->SetSelectedComponent(0);
	threshold->ThresholdBetween(_index, _index);

	vtkDataSetSurfaceFilter* surface = vtkDataSetSurfaceFilter::New();
	surface->SetInputConnection(threshold->GetOutputPort());

	if (_type == GEOLIB::POINT)
	{
		VtkCompositePointToGlyphFilter* pts = new VtkCompositePointToGlyphFilter(surface);
	}
	else if (_type == GEOLIB::POLYLINE)
	{
		VtkCompositeLineToTubeFilter* ltt = new VtkCompositeLineToTubeFilter(surface);
	}

	_outputAlgorithm = surface;
}

