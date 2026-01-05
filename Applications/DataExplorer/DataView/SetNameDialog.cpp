// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "SetNameDialog.h"

#include <QDialogButtonBox>
#include <QLabel>
#include <QLineEdit>
#include <QVBoxLayout>

SetNameDialog::SetNameDialog(const std::string& geo_object_type, std::size_t id,
                             const std::string& old_name = "", QDialog* parent)
    : QDialog(parent)
{
    QString const& label =
        QString::fromStdString(geo_object_type) + "#" + QString::number(id);
    setupDialog(label, old_name);
    show();
}

SetNameDialog::~SetNameDialog()
{
    delete _buttonBox;
    delete _layout;
    delete _new_name;
    delete _txt_label;
}

void SetNameDialog::setupDialog(const QString& label,
                                const std::string& old_name)
{
    _layout = new QVBoxLayout(this);
    QString dialog_text("Please enter a name for " + label);
    _txt_label = new QLabel(dialog_text, this);
    _new_name = new QLineEdit(QString::fromStdString(old_name));

    setWindowTitle("Set name...");
    _layout->addWidget(_txt_label);
    _layout->addWidget(_new_name);
    _buttonBox =
        new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    connect(_buttonBox, SIGNAL(accepted()), this, SLOT(accept()));
    connect(_buttonBox, SIGNAL(rejected()), this, SLOT(reject()));
    _layout->addWidget(_buttonBox);

    setLayout(_layout);
}

std::string SetNameDialog::getNewName()
{
    return _new_name->text().toStdString();
}

void SetNameDialog::accept()
{
    this->done(QDialog::Accepted);
}

void SetNameDialog::reject()
{
    this->done(QDialog::Rejected);
}
