/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file QVtkDataSetMapper.cpp
 *
 * Created on 2010-11-12 by Lars Bilke
 */

// ** INCLUDES **
#include "QVtkDataSetMapper.h"

#include <vtkObjectFactory.h>

vtkStandardNewMacro(QVtkDataSetMapper);

QVtkDataSetMapper::QVtkDataSetMapper()
	: QObject(NULL)
{
}

QVtkDataSetMapper::~QVtkDataSetMapper()
{
}

void QVtkDataSetMapper::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
}

void QVtkDataSetMapper::SetScalarVisibility( bool on )
{
	vtkDataSetMapper::SetScalarVisibility((int)on);
}

