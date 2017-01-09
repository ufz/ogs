/**
 * \file
 * \author Lars Bilke
 * \date   2010-11-12
 * \brief  Implementation of the QVtkDataSetMapper class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
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

