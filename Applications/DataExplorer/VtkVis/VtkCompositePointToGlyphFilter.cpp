/**
 * \file
 * \author Lars Bilke
 * \date   2010-10-21
 * \brief  Implementation of the VtkCompositePointToGlyphFilter class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkCompositePointToGlyphFilter.h"

#include <vtkDataSetAlgorithm.h>
#include <vtkGlyph3D.h>
#include <vtkPointData.h>
#include <vtkSphereSource.h>

VtkCompositePointToGlyphFilter::VtkCompositePointToGlyphFilter( vtkAlgorithm* inputAlgorithm )
    : VtkCompositeFilter(inputAlgorithm), glyphSource_(nullptr)
{
    this->init();
}

VtkCompositePointToGlyphFilter::~VtkCompositePointToGlyphFilter()
{
    glyphSource_->Delete();
}

void VtkCompositePointToGlyphFilter::init()
{
    this->inputDataObjectType_ = VTK_DATA_SET;
    this->outputDataObjectType_ = VTK_POLY_DATA;

    std::size_t nPoints = static_cast<vtkDataSetAlgorithm*>(inputAlgorithm_)
                        ->GetOutput()->GetPointData()->GetNumberOfTuples();
    int phi (10 - static_cast<std::size_t>(nPoints / 2000.0));
    int theta (phi);
    if (phi < 4)
    {
        phi = 4;
        theta = 4; // for theta 3 would be possible, too, but 4 looks much better
    }

    double default_radius(GetInitialRadius());
    glyphSource_ = vtkSphereSource::New();
    glyphSource_->SetRadius(default_radius);
    glyphSource_->SetPhiResolution(phi);
    glyphSource_->SetThetaResolution(theta);
    (*algorithmUserProperties_)["Radius"] = default_radius;

    (*algorithmUserProperties_)["PhiResolution"] = phi;
    (*algorithmUserProperties_)["ThetaResolution"] = theta;

    vtkGlyph3D* glyphFilter = vtkGlyph3D::New();
    glyphFilter->ScalingOn();   // KR important to scale glyphs with double precision (e.g. 0.1 of their size for small datasets)
    //glyphFilter->SetScaleModeToScaleByScalar();  // KR can easily obscure view when scalar values have large differences (this is also the default scaling method)
    glyphFilter->SetScaleModeToDataScalingOff(); // KR scaling is possible but scalar values are ignored
    glyphFilter->SetScaleFactor(1.0);
    glyphFilter->SetSourceConnection(glyphSource_->GetOutputPort());
    glyphFilter->SetInputConnection(inputAlgorithm_->GetOutputPort());
    //(*algorithmUserProperties_)["ScaleMode"] = 0;
    //(*algorithmUserProperties_)["ScaleFactor"] = 1.0;
    //(*algorithmUserProperties_)["ColorMode"] = glyphFilter->GetColorMode();
    //(*algorithmUserProperties_)["VectorMode"] = glyphFilter->GetVectorMode();
    //(*algorithmUserProperties_)["Orient"] = glyphFilter->GetOrient();

    outputAlgorithm_ = glyphFilter;
}

void VtkCompositePointToGlyphFilter::SetUserProperty( QString name, QVariant value )
{
    VtkAlgorithmProperties::SetUserProperty(name, value);

    if (name.compare("Radius") == 0)
    {
        glyphSource_->SetRadius(value.toDouble());
    }
    else if (name.compare("PhiResolution") == 0)
    {
        glyphSource_->SetPhiResolution(value.toInt());
    }
    else if (name.compare("ThetaResolution") == 0)
    {
        glyphSource_->SetThetaResolution(value.toInt());
    }
    else if (name.compare("ScaleMode") == 0)
    {
        static_cast<vtkGlyph3D*>(outputAlgorithm_)->SetScaleMode(value.toInt());
    }
    else if (name.compare("ScaleFactor") == 0)
    {
        static_cast<vtkGlyph3D*>(outputAlgorithm_)->SetScaleFactor(value.toDouble());
    }
    else if (name.compare("ColorMode") == 0)
    {
        static_cast<vtkGlyph3D*>(outputAlgorithm_)->SetColorMode(value.toInt());
    }
    else if (name.compare("VectorMode") == 0)
    {
        static_cast<vtkGlyph3D*>(outputAlgorithm_)->SetVectorMode(value.toInt());
    }
    else if (name.compare("Orient") == 0)
    {
        static_cast<vtkGlyph3D*>(outputAlgorithm_)->SetOrient(value.toBool());
    }
}
