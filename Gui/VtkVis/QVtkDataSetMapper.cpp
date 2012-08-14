/**
 * \file QVtkDataSetMapper.cpp
 * 12/11/2010 LB Initial implementation
 *
 * Implementation of QVtkDataSetMapper class
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

