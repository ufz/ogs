/**
 * \file
 * \author Karsten Rink
 * \date   2010-03-19
 * \brief  Implementation of the VtkMeshSource class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "VtkMeshSource.h"

#include "logog/include/logog.hpp"

#include "Color.h"

#include <vtkProperty.h>


vtkStandardNewMacro(VtkMeshSource);

VtkMeshSource::VtkMeshSource() :
		_grid(nullptr), _matName("MaterialIDs")
{
	_removable = false; // From VtkAlgorithmProperties
	this->SetNumberOfInputPorts(0);

	const GeoLib::Color* c = GeoLib::getRandomColor();
	vtkProperty* vtkProps = GetProperties();
	vtkProps->SetColor((*c)[0] / 255.0,(*c)[1] / 255.0,(*c)[2] / 255.0);
	delete c;
	vtkProps->SetEdgeVisibility(1);
}

VtkMeshSource::~VtkMeshSource()
{
}

void VtkMeshSource::PrintSelf( ostream& os, vtkIndent indent )
{
	this->Superclass::PrintSelf(os,indent);
}

int VtkMeshSource::RequestData( vtkInformation* request,
                                vtkInformationVector** inputVector,
                                vtkInformationVector* outputVector )
{
	return this->Superclass::RequestData(request, inputVector, outputVector);
}

void VtkMeshSource::SetUserProperty( QString name, QVariant value )
{
	VtkAlgorithmProperties::SetUserProperty(name, value);

	(*_algorithmUserProperties)[name] = value;
}
