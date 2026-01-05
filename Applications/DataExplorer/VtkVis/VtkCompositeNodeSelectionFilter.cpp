// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "VtkCompositeNodeSelectionFilter.h"

#include <vtkDataSetAlgorithm.h>
#include <vtkGlyph3D.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>

#include <memory>

#include "VtkPointsSource.h"

VtkCompositeNodeSelectionFilter::VtkCompositeNodeSelectionFilter(
    vtkAlgorithm* inputAlgorithm)
    : VtkCompositeFilter(inputAlgorithm)
{
    // this->init();
}

VtkCompositeNodeSelectionFilter::~VtkCompositeNodeSelectionFilter()
{
    for (auto& item : _selection)
    {
        delete item;
    }
}

void VtkCompositeNodeSelectionFilter::init()
{
    this->_inputDataObjectType = VTK_DATA_SET;
    this->_outputDataObjectType = VTK_POLY_DATA;

    if (!_selection.empty())
    {
        vtkSmartPointer<VtkPointsSource> point_source =
            vtkSmartPointer<VtkPointsSource>::New();
        point_source->setPoints(&_selection);

        vtkSmartPointer<vtkSphereSource> _glyphSource =
            vtkSmartPointer<vtkSphereSource>::New();
        _glyphSource->SetRadius(this->GetInitialRadius());

        vtkGlyph3D* glyphFilter = vtkGlyph3D::New();
        glyphFilter->SetSourceConnection(_glyphSource->GetOutputPort());
        glyphFilter->SetInputConnection(point_source->GetOutputPort());

        _outputAlgorithm = glyphFilter;
    }
    else
    {
        _outputAlgorithm = nullptr;
    }
}

void VtkCompositeNodeSelectionFilter::setSelectionArray(
    const std::vector<unsigned>& point_indeces)
{
    for (unsigned int point_index : point_indeces)
    {
        double* coords = static_cast<vtkDataSetAlgorithm*>(_inputAlgorithm)
                             ->GetOutput()
                             ->GetPoint(point_index);
        auto* p(new GeoLib::Point(coords[0], coords[1], coords[2]));
        _selection.push_back(p);
    }
    init();
}
