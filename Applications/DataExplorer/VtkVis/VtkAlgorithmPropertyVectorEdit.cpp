// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

// ** INCLUDES **
#include "VtkAlgorithmPropertyVectorEdit.h"

#include <QDoubleValidator>
#include <QHBoxLayout>
#include <QIntValidator>
#include <QLineEdit>
#include <QSize>
#include <utility>

#include "VtkAlgorithmProperties.h"

VtkAlgorithmPropertyVectorEdit::VtkAlgorithmPropertyVectorEdit(
    const QList<QString> contents,
    QString name,
    QVariant::Type type,
    VtkAlgorithmProperties* algProps,
    QWidget* parent /*= 0*/)
    : QWidget(parent), _name(std::move(name)), _algProps(algProps), _type(type)
{
    auto* layout = new QHBoxLayout;
    layout->setSpacing(3);
    layout->setContentsMargins(0, 0, 0, 0);

    foreach (QString content, contents)
    {
        auto* lineEdit = new QLineEdit(content, this);
        layout->addWidget(lineEdit);

        switch (_type)
        {
            case QVariant::Double:
                lineEdit->setValidator(new QDoubleValidator(this));
                break;

            case QVariant::Int:
                lineEdit->setValidator(new QIntValidator(this));
                break;

            default:
                break;
        }

        connect(lineEdit, SIGNAL(editingFinished()), this, SLOT(setNewValue()));
    }

    this->setLayout(layout);
    this->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);
}

VtkAlgorithmPropertyVectorEdit::~VtkAlgorithmPropertyVectorEdit() = default;

void VtkAlgorithmPropertyVectorEdit::setNewValue()
{
    QLayout* layout = this->layout();
    QList<QVariant> list;
    for (int i = 0; i < layout->count(); ++i)
    {
        auto* lineEdit = static_cast<QLineEdit*>(layout->itemAt(i)->widget());
        list.push_back(QVariant(lineEdit->text()));
    }

    _algProps->SetUserVectorProperty(_name, list);

    emit editingFinished();
}
