/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file VtkAlgorithmPropertyCheckbox.cpp
 *
 * Created on 2010-10-20 by Lars Bilke
 */

// ** INCLUDES **
#include "VtkAlgorithmPropertyCheckbox.h"

#include "VtkAlgorithmProperties.h"

VtkAlgorithmPropertyCheckbox::VtkAlgorithmPropertyCheckbox(const bool value,
                                                           const QString& name,
                                                           VtkAlgorithmProperties* algProps,
                                                           QWidget* parent /*= 0*/ )
	: QCheckBox(parent), _name(name), _algProps(algProps)
{
	this->setChecked(value);
	connect(this, SIGNAL(stateChanged(int)), this, SLOT(setNewValue(int)));
}

VtkAlgorithmPropertyCheckbox::~VtkAlgorithmPropertyCheckbox()
{
}

void VtkAlgorithmPropertyCheckbox::setNewValue( int state )
{
	bool boolState = (bool)state;
	_algProps->SetUserProperty(_name, QVariant(boolState));
}
