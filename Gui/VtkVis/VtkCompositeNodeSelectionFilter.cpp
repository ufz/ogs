/**
 * \file
 * \author Karsten Rink
 * \date   2013-04-16
 * \brief  Implementation of the VtkCompositeNodeSelectionFilter class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkCompositeNodeSelectionFilter.h"
#include "VtkCompositePointToGlyphFilter.h"
#include "VtkPointsSource.h"

#include <vtkDataSetAlgorithm.h>
#include <vtkSmartPointer.h>


VtkCompositeNodeSelectionFilter::VtkCompositeNodeSelectionFilter( vtkAlgorithm* inputAlgorithm )
	: VtkCompositeFilter(inputAlgorithm)
{
	int a = 3;
	//this->init();
}

VtkCompositeNodeSelectionFilter::~VtkCompositeNodeSelectionFilter()
{
	for (unsigned i=0; i<_selection.size(); ++i)
		delete _selection[i];
}

void VtkCompositeNodeSelectionFilter::init()
{
	this->_inputDataObjectType = VTK_DATA_SET;
	this->_outputDataObjectType = VTK_POLY_DATA;

	if (!_selection.empty())
	{
		VtkPointsSource* point_source = VtkPointsSource::New();
		point_source->setPoints(&_selection);

		VtkCompositeFilter* glyphs = new VtkCompositePointToGlyphFilter(point_source);
		glyphs->SetUserProperty("Radius", this->GetInitialRadius());
		_outputAlgorithm = glyphs->GetOutputAlgorithm();
	}
	else
		_outputAlgorithm = nullptr;
}

void VtkCompositeNodeSelectionFilter::setSelectionArray(const std::vector<unsigned> &point_indeces)
{
	for (unsigned i=0; i<point_indeces.size(); ++i)
	{
		unsigned idx = point_indeces[i];
		double * coords = static_cast<vtkDataSetAlgorithm*>(_inputAlgorithm)->GetOutput()->GetPoint(point_indeces[i]);
		GeoLib::Point* p (new GeoLib::Point(coords));
		_selection.push_back(p);
	}
	init(); 
}

