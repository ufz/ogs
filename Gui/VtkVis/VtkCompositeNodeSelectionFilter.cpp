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
#include <memory>

#include "VtkCompositeNodeSelectionFilter.h"
//#include "VtkCompositePointToGlyphFilter.h"
#include "VtkPointsSource.h"

#include <vtkDataSetAlgorithm.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkGlyph3D.h>


VtkCompositeNodeSelectionFilter::VtkCompositeNodeSelectionFilter( vtkAlgorithm* inputAlgorithm )
: VtkCompositeFilter(inputAlgorithm)
{
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
		vtkSmartPointer<VtkPointsSource> point_source = vtkSmartPointer<VtkPointsSource>::New();
		point_source->setPoints(&_selection);

		vtkSmartPointer<vtkSphereSource> _glyphSource = vtkSmartPointer<vtkSphereSource>::New();
			_glyphSource->SetRadius(this->GetInitialRadius());

		vtkGlyph3D* glyphFilter = vtkGlyph3D::New();
			glyphFilter->SetSourceData(_glyphSource->GetOutput());
			glyphFilter->SetInputConnection(point_source->GetOutputPort());

		_outputAlgorithm = glyphFilter;
	}
	else
		_outputAlgorithm = nullptr;
}

void VtkCompositeNodeSelectionFilter::setSelectionArray(const std::vector<unsigned> &point_indeces)
{
	for (unsigned i=0; i<point_indeces.size(); ++i)
	{
		double * coords = static_cast<vtkDataSetAlgorithm*>(_inputAlgorithm)->GetOutput()->GetPoint(point_indeces[i]);
		GeoLib::Point* p (new GeoLib::Point(coords));
		_selection.push_back(p);
	}
	init(); 
}

