/**
 * \file VtkCompositeFilter.cpp
 * 19/10/2010 LB Initial implementation
 * 
 * Implementation of VtkCompositeFilter class
 */

// ** INCLUDES **
#include "VtkCompositeFilter.h"

#include "VtkCompositeImageToCylindersFilter.h"

#include <QVector>
#include <QString>
#include <QMapIterator>

VtkCompositeFilter::VtkCompositeFilter(vtkAlgorithm* inputAlgorithm)
: _inputAlgorithm(inputAlgorithm)
{
}

VtkCompositeFilter::~VtkCompositeFilter()
{
	 _outputAlgorithm->Delete();
}
