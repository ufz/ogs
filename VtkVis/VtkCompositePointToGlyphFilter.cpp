/**
 * \file VtkCompositePointToGlyphFilter.cpp
 * 21/10/2010 LB Initial implementation
 * 
 * Implementation of VtkCompositePointToGlyphFilter class
 */

// ** INCLUDES **
#include "VtkCompositePointToGlyphFilter.h"

#include <vtkSphereSource.h>
#include <vtkGlyph3D.h>

VtkCompositePointToGlyphFilter::VtkCompositePointToGlyphFilter( vtkAlgorithm* inputAlgorithm )
: VtkCompositeFilter(inputAlgorithm)
{
	this->init();
}

VtkCompositePointToGlyphFilter::~VtkCompositePointToGlyphFilter()
{
	_glyphSource->Delete();
}

void VtkCompositePointToGlyphFilter::init()
{
	this->_inputDataObjectType = VTK_DATA_SET;
	this->_outputDataObjectType = VTK_POLY_DATA;

	int default_radius(150);
	_glyphSource = vtkSphereSource::New();
	_glyphSource->SetRadius(default_radius);
	_glyphSource->SetPhiResolution(10);
	_glyphSource->SetThetaResolution(10);
	(*_algorithmUserProperties)["Radius"] = default_radius;
	(*_algorithmUserProperties)["PhiResolution"] = 10;
	(*_algorithmUserProperties)["ThetaResolution"] = 10;

	vtkGlyph3D* glyphFilter = vtkGlyph3D::New();
	glyphFilter->ScalingOn();		// KR important to scale glyphs with double precision (e.g. 0.1 of their size for small datasets)
	//glyphFilter->SetScaleModeToScaleByScalar();  // KR can easily obscure view when scalar values have large differences (this is also the default scaling method)
	glyphFilter->SetScaleModeToDataScalingOff();	// KR scaling is possible but scalar values are ignored
	glyphFilter->SetScaleFactor(1.0);
	glyphFilter->SetSource(_glyphSource->GetOutput());
	glyphFilter->SetInputConnection(_inputAlgorithm->GetOutputPort());
	//(*_algorithmUserProperties)["ScaleMode"] = 0;
	(*_algorithmUserProperties)["ScaleFactor"] = 1.0;
	//(*_algorithmUserProperties)["ColorMode"] = glyphFilter->GetColorMode();
	//(*_algorithmUserProperties)["VectorMode"] = glyphFilter->GetVectorMode();
	//(*_algorithmUserProperties)["Orient"] = glyphFilter->GetOrient();

	_outputAlgorithm = glyphFilter;
}

void VtkCompositePointToGlyphFilter::SetUserProperty( QString name, QVariant value )
{
	VtkAlgorithmProperties::SetUserProperty(name, value);

	if (name.compare("Radius") == 0)
		_glyphSource->SetRadius(value.toDouble());
	else if (name.compare("PhiResolution") == 0)
		_glyphSource->SetPhiResolution(value.toInt());
	else if (name.compare("ThetaResolution") == 0)
		_glyphSource->SetThetaResolution(value.toInt());
	else if (name.compare("ScaleMode") == 0)
		static_cast<vtkGlyph3D*>(_outputAlgorithm)->SetScaleMode(value.toInt());
	else if (name.compare("ScaleFactor") == 0)
		static_cast<vtkGlyph3D*>(_outputAlgorithm)->SetScaleFactor(value.toDouble());
	else if (name.compare("ColorMode") == 0)
		static_cast<vtkGlyph3D*>(_outputAlgorithm)->SetColorMode(value.toInt());
	else if (name.compare("VectorMode") == 0)
		static_cast<vtkGlyph3D*>(_outputAlgorithm)->SetVectorMode(value.toInt());
	else if (name.compare("Orient") == 0)
		static_cast<vtkGlyph3D*>(_outputAlgorithm)->SetOrient(value.toBool());
}
