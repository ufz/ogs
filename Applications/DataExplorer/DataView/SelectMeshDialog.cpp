/**
 * \file
 * \author Karsten Rink
 * \date   2012-04-20
 * \brief  Implementation of the SelectMeshDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
    QDialog(parent), _geo_object(geo_object)
{
    setupDialog(msh_names);
    show();
}

SelectMeshDialog::~SelectMeshDialog()
{
    delete _buttonBox;
    delete _layout;
    delete _msh_names;
    delete _txt_label;
}

void SelectMeshDialog::setupDialog(const std::list<std::string> &msh_names)
{
    _layout = new QVBoxLayout(this);
    QString dialog_text("Select Mesh");
    _txt_label = new QLabel(this);
    _txt_label->setText(dialog_text);


    _msh_names = new QComboBox();
    for (auto it = msh_names.begin(); it != msh_names.end(); ++it)
        _msh_names->addItem(QString::fromStdString(*it));

    setWindowTitle("Select Mesh...");
    _layout->addWidget( _txt_label );
    _layout->addWidget( _msh_names );
    _buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    connect(_buttonBox, SIGNAL(accepted()), this, SLOT(accept()));
    connect(_buttonBox, SIGNAL(rejected()), this, SLOT(reject()));
    _layout->addWidget( _buttonBox );

    setLayout(_layout);
}

void SelectMeshDialog::accept()
{
    //emit requestNameChange(_parent_name, GeoLib::convertGeoType(_object_type_name), _id, _new_name->text().toStdString());
    this->done(QDialog::Accepted);
}

void SelectMeshDialog::reject()
{
    this->done(QDialog::Rejected);
}
