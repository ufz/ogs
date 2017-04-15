/**
 * \file
 * \author Lars Bilke
 * \date   2010-10-20
 * \brief  Implementation of the VtkAlgorithmPropertyCheckbox class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
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

VtkAlgorithmPropertyCheckbox::~VtkAlgorithmPropertyCheckbox() = default;

void VtkAlgorithmPropertyCheckbox::setNewValue( int state )
{
    auto boolState = (bool)state;
    _algProps->SetUserProperty(_name, QVariant(boolState));
}
