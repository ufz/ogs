/**
 * \file
 * \author Karsten Rink
 * \date   2011-10-26
 * \brief  Implementation of the SetNameDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SetNameDialog.h"

#include <QDialogButtonBox>
#include <QDialogButtonBox>
#include <QLabel>
#include <QLineEdit>
#include <QVBoxLayout>

SetNameDialog::SetNameDialog(const std::string &geo_object_type, std::size_t id, const std::string &old_name = "", QDialog* parent)
:QDialog(parent)
{
    QString const& label = QString::fromStdString(geo_object_type) + "#" + QString::number(id);
    setupDialog(label, old_name);
    show();
}

SetNameDialog::~SetNameDialog()
{
    delete buttonBox_;
    delete layout_;
    delete new_name_;
    delete txt_label_;
}

void SetNameDialog::setupDialog(const QString &label, const std::string &old_name)
{
    layout_ = new QVBoxLayout(this);
    QString dialog_text("Please enter a name for " + label);
    txt_label_ = new QLabel(dialog_text, this);
    new_name_ = new QLineEdit(QString::fromStdString(old_name));

    setWindowTitle("Set name...");
    layout_->addWidget( txt_label_ );
    layout_->addWidget( new_name_ );
    buttonBox_ = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    connect(buttonBox_, SIGNAL(accepted()), this, SLOT(accept()));
    connect(buttonBox_, SIGNAL(rejected()), this, SLOT(reject()));
    layout_->addWidget( buttonBox_ );

    setLayout(layout_);
}

std::string SetNameDialog::getNewName()
{
    return new_name_->text().toStdString();
}

void SetNameDialog::accept()
{
    this->done(QDialog::Accepted);
}

void SetNameDialog::reject()
{
    this->done(QDialog::Rejected);
}
