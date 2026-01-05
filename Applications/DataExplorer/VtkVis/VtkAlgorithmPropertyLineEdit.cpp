// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

// ** INCLUDES **
#include "VtkAlgorithmPropertyLineEdit.h"

#include <QDoubleValidator>
#include <QIntValidator>
#include <utility>

#include "VtkAlgorithmProperties.h"

VtkAlgorithmPropertyLineEdit::VtkAlgorithmPropertyLineEdit(
    const QString& contents,
    QString name,
    QVariant::Type type,
    VtkAlgorithmProperties* algProps,
    QWidget* parent /*= 0*/)
    : QLineEdit(contents, parent),
      _name(std::move(name)),
      _algProps(algProps),
      _type(type)
{
    switch (_type)
    {
        case QVariant::Double:
            this->setValidator(new QDoubleValidator(this));
            break;

        case QVariant::Int:
            this->setValidator(new QIntValidator(this));
            break;

        default:
            break;
    }

    connect(this, SIGNAL(editingFinished()), this, SLOT(setNewValue()));
}

VtkAlgorithmPropertyLineEdit::~VtkAlgorithmPropertyLineEdit() = default;

void VtkAlgorithmPropertyLineEdit::setNewValue()
{
    QVariant value(this->text());
    if (value.convert(_type))
    {
        _algProps->SetUserProperty(_name, value);
    }
}
