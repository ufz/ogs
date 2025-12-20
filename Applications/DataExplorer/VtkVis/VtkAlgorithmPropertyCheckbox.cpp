// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

// ** INCLUDES **
#include "VtkAlgorithmPropertyCheckbox.h"

#include <utility>

#include "VtkAlgorithmProperties.h"

VtkAlgorithmPropertyCheckbox::VtkAlgorithmPropertyCheckbox(
    const bool value,
    QString name,
    VtkAlgorithmProperties* algProps,
    QWidget* parent /*= 0*/)
    : QCheckBox(parent), _name(std::move(name)), _algProps(algProps)
{
    this->setChecked(value);
    connect(this, SIGNAL(stateChanged(int)), this, SLOT(setNewValue(int)));
}

VtkAlgorithmPropertyCheckbox::~VtkAlgorithmPropertyCheckbox() = default;

void VtkAlgorithmPropertyCheckbox::setNewValue(int state)
{
    auto boolState = static_cast<bool>(state);
    _algProps->SetUserProperty(_name, QVariant(boolState));
}
