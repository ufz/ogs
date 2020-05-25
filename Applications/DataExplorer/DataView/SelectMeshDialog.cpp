/**
 * \file
 * \author Karsten Rink
 * \date   2012-04-20
 * \brief  Implementation of the SelectMeshDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SelectMeshDialog.h"
#include "GeoObject.h"

#include <QDialogButtonBox>
#include <QLabel>
#include <QComboBox>
#include <QVBoxLayout>

SelectMeshDialog::SelectMeshDialog(const GeoLib::GeoObject* geo_object, const std::list<std::string> &msh_names, QDialog* parent) :
    QDialog(parent), geo_object_(geo_object)
{
    setupDialog(msh_names);
    show();
}

SelectMeshDialog::~SelectMeshDialog()
{
    delete buttonBox_;
    delete layout_;
    delete msh_names_;
    delete txt_label_;
}

void SelectMeshDialog::setupDialog(const std::list<std::string> &msh_names)
{
    layout_ = new QVBoxLayout(this);
    QString dialog_text("Select Mesh");
    txt_label_ = new QLabel(this);
    txt_label_->setText(dialog_text);


    msh_names_ = new QComboBox();
    for (const auto& msh_name : msh_names)
    {
        msh_names_->addItem(QString::fromStdString(msh_name));
    }

    setWindowTitle("Select Mesh...");
    layout_->addWidget( txt_label_ );
    layout_->addWidget( msh_names_ );
    buttonBox_ = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    connect(buttonBox_, SIGNAL(accepted()), this, SLOT(accept()));
    connect(buttonBox_, SIGNAL(rejected()), this, SLOT(reject()));
    layout_->addWidget( buttonBox_ );

    setLayout(layout_);
}

void SelectMeshDialog::accept()
{
    //emit requestNameChange(parent_name_, GeoLib::convertGeoType(object_type_name_), id_, new_name_->text().toStdString());
    this->done(QDialog::Accepted);
}

void SelectMeshDialog::reject()
{
    this->done(QDialog::Rejected);
}
