/**
 * \file VtkCompositeFilter.cpp
 * 19/10/2010 LB Initial implementation
 *
 * Implementation of VtkCompositeFilter class
 */

// ** INCLUDES **
#include "VtkCompositeFilter.h"

#include <vtkAlgorithm.h>

#include <QMapIterator>
#include <QString>
#include <QVector>

VtkCompositeFilter::VtkCompositeFilter(vtkAlgorithm* inputAlgorithm)
	: _inputDataObjectType(0), _outputDataObjectType(1),
	  _inputAlgorithm(inputAlgorithm)
{
}

VtkCompositeFilter::~VtkCompositeFilter()
{
	_outputAlgorithm->Delete();
}
