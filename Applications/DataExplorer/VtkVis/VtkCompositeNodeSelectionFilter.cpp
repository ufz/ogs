/**
 * \file
 * \author Karsten Rink
 * \date   2013-04-16
 * \brief  Implementation of the VtkCompositeNodeSelectionFilter class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
    for (auto& item : selection_)
    {
        delete item;
    }
}

void VtkCompositeNodeSelectionFilter::init()
{
    this->inputDataObjectType_ = VTK_DATA_SET;
    this->outputDataObjectType_ = VTK_POLY_DATA;

    if (!selection_.empty())
    {
        vtkSmartPointer<VtkPointsSource> point_source = vtkSmartPointer<VtkPointsSource>::New();
        point_source->setPoints(&selection_);

        vtkSmartPointer<vtkSphereSource> glyphSource_ = vtkSmartPointer<vtkSphereSource>::New();
            glyphSource_->SetRadius(this->GetInitialRadius());

        vtkGlyph3D* glyphFilter = vtkGlyph3D::New();
            glyphFilter->SetSourceConnection(glyphSource_->GetOutputPort());
            glyphFilter->SetInputConnection(point_source->GetOutputPort());

        outputAlgorithm_ = glyphFilter;
    }
    else
    {
        outputAlgorithm_ = nullptr;
    }
}

void VtkCompositeNodeSelectionFilter::setSelectionArray(const std::vector<unsigned> &point_indeces)
{
    for (unsigned int point_index : point_indeces)
    {
        double* coords = static_cast<vtkDataSetAlgorithm*>(inputAlgorithm_)
                             ->GetOutput()
                             ->GetPoint(point_index);
        auto* p(new GeoLib::Point(coords[0], coords[1], coords[2]));
        selection_.push_back(p);
    }
    init();
}

