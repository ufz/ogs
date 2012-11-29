/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file VtkCompositeFilter.cpp
 *
 * Created on 2010-10-19 by Lars Bilke
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
