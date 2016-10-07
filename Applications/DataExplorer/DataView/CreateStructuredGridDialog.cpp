/**
 * \file   CreateStructuredGridDialog.cpp
 * \author Karsten Rink
 * \date   2016-02-04
 * \brief  Implementation of the CreateStructuredGridDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateStructuredGridDialog.h"

#include <QIntValidator>

#include "GeoLib/Point.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"

#include "StrictDoubleValidator.h"
#include "OGSError.h"


CreateStructuredGridDialog::CreateStructuredGridDialog(QDialog* parent) : QDialog(parent)
{
    setupUi(this);
    setValidators();
}

void CreateStructuredGridDialog::setValidators()
{
    StrictDoubleValidator* origin_x_validator = new StrictDoubleValidator(this);
    this->xOriginEdit->setValidator (origin_x_validator);
    StrictDoubleValidator* origin_y_validator = new StrictDoubleValidator(this);
    this->yOriginEdit->setValidator (origin_y_validator);
    StrictDoubleValidator* origin_z_validator = new StrictDoubleValidator(this);
    this->zOriginEdit->setValidator (origin_z_validator);

    StrictDoubleValidator* x_length_validator = new StrictDoubleValidator(0, 10000000, 10, this);
    this->xLengthEdit->setValidator (x_length_validator);
    StrictDoubleValidator* y_length_validator = new StrictDoubleValidator(0, 10000000, 10, this);
    this->yLengthEdit->setValidator (y_length_validator);
    StrictDoubleValidator* z_length_validator = new StrictDoubleValidator(0, 10000000, 10, this);
    this->zLengthEdit->setValidator (z_length_validator);

    QIntValidator* x_n_elem_validator = new QIntValidator(1, 10000000, this);
    this->xElemEdit->setValidator (x_n_elem_validator);
    QIntValidator* y_n_elem_validator = new QIntValidator(1, 10000000, this);
    this->yElemEdit->setValidator (y_n_elem_validator);
    QIntValidator* z_n_elem_validator = new QIntValidator(1, 10000000, this);
    this->zElemEdit->setValidator (z_n_elem_validator);
}

void CreateStructuredGridDialog::on_lineButton_toggled() const
{
    this->yLengthLabel->setEnabled(false);
    this->yLengthEdit->setEnabled(false);
    this->zLengthLabel->setEnabled(false);
    this->zLengthEdit->setEnabled(false);
    this->yElemLabel->setEnabled(false);
    this->yElemEdit->setEnabled(false);
    this->zElemLabel->setEnabled(false);
    this->zElemEdit->setEnabled(false);
}

void CreateStructuredGridDialog::enable2dWidgets() const
{
    this->yLengthLabel->setEnabled(true);
    this->yLengthEdit->setEnabled(true);
    this->zLengthLabel->setEnabled(false);
    this->zLengthEdit->setEnabled(false);
    this->yElemLabel->setEnabled(true);
    this->yElemEdit->setEnabled(true);
    this->zElemLabel->setEnabled(false);
    this->zElemEdit->setEnabled(false);
}

void CreateStructuredGridDialog::enable3dWidgets() const
{
    this->yLengthLabel->setEnabled(true);
    this->yLengthEdit->setEnabled(true);
    this->zLengthLabel->setEnabled(true);
    this->zLengthEdit->setEnabled(true);
    this->yElemLabel->setEnabled(true);
    this->yElemEdit->setEnabled(true);
    this->zElemLabel->setEnabled(true);
    this->zElemEdit->setEnabled(true);
}

void CreateStructuredGridDialog::on_meshExtentButton_toggled()
{
    this->xLengthLabel->setText("Mesh size in x");
    this->yLengthLabel->setText("Mesh size in y");
    this->zLengthLabel->setText("Mesh size in z");
}

void CreateStructuredGridDialog::on_elemExtentButton_toggled()
{
    this->xLengthLabel->setText("Element size in x");
    this->yLengthLabel->setText("Element size in y");
    this->zLengthLabel->setText("Element size in z");
}

bool CreateStructuredGridDialog::inputIsEmpty() const
{
    QString const type_str = (this->meshExtentButton->isChecked()) ? "mesh" : "element";
    if (this->xLengthEdit->text().isEmpty())
    {
        OGSError::box("Please specify " + type_str + "\nextent in x-direction.");
        return true;
    }
    if (this->xElemEdit->text().isEmpty())
    {
        OGSError::box("Please specify number of\nelements in x-direction.");
        return true;
    }
    if (this->xOriginEdit->text().isEmpty() ||
        this->yOriginEdit->text().isEmpty() ||
        this->zOriginEdit->text().isEmpty())
    {
        OGSError::box("Please specify coordinates\nof mesh origin.");
        return true;
    }
    if (this->meshNameEdit->text().isEmpty())
    {
        OGSError::box("Please specify mesh name.");
        return true;
    }

    if (!this->lineButton->isChecked())
    {
        if (this->yLengthEdit->text().isEmpty())
        {
            OGSError::box("Please specify " + type_str + "\nextent in y-direction.");
            return true;
        }
        if (this->yElemEdit->text().isEmpty())
        {
            OGSError::box("Please specify number of\nelements in y-direction.");
            return true;
        }
    }

    if (this->prismButton->isChecked() || this->hexButton->isChecked())
    {
        if (this->zLengthEdit->text().isEmpty())
        {
            OGSError::box("Please specify " + type_str + "\nextent in z-direction.");
            return true;
        }
        if (this->zElemEdit->text().isEmpty())
        {
            OGSError::box("Please specify number of\nelements in z-direction.");
            return true;
        }
    }
    return false;
}

void CreateStructuredGridDialog::accept()
{
    if (inputIsEmpty())
        return;

    if ((this->xLengthEdit->text().toDouble() <= 0) ||
        (this->yLengthEdit->text().toDouble() <= 0) ||
        (this->zLengthEdit->text().toDouble() <= 0))
    {
        OGSError::box("Length needs to be larger than 0.");
        return;
    }

    if ((this->xElemEdit->text().toDouble() <= 0) ||
        (this->yElemEdit->text().toDouble() <= 0) ||
        (this->zElemEdit->text().toDouble() <= 0))
    {
        OGSError::box("Number of elements needs to be larger than 0.");
        return;
    }

    GeoLib::Point const origin(this->xOriginEdit->text().toDouble(),
                               this->yOriginEdit->text().toDouble(),
                               this->zOriginEdit->text().toDouble());
    std::string const name (this->meshNameEdit->text().toStdString());
    MeshLib::Mesh* mesh (nullptr);
    if (this->lineButton->isChecked())
        if (this->meshExtentButton->isChecked())
            mesh = MeshLib::MeshGenerator::generateLineMesh(
                this->xLengthEdit->text().toDouble(), this->xElemEdit->text().toInt(), origin, name);
        else
            mesh = MeshLib::MeshGenerator::generateLineMesh(
                this->xElemEdit->text().toInt(), this->xLengthEdit->text().toDouble(), origin, name);
    else if (this->triButton->isChecked())
        if (this->meshExtentButton->isChecked())
            mesh = MeshLib::MeshGenerator::generateRegularTriMesh(
                this->xLengthEdit->text().toDouble(), this->yLengthEdit->text().toDouble(),
                this->xElemEdit->text().toInt(), this->yElemEdit->text().toInt(),
                origin, name);
        else
            mesh = MeshLib::MeshGenerator::generateRegularTriMesh(
                this->xElemEdit->text().toInt(), this->yElemEdit->text().toInt(),
                this->xLengthEdit->text().toDouble(), this->yLengthEdit->text().toDouble(),
                origin, name);
    else if (this->quadButton->isChecked())
        if (this->meshExtentButton->isChecked())
            mesh = MeshLib::MeshGenerator::generateRegularQuadMesh(
                this->xLengthEdit->text().toDouble(), this->yLengthEdit->text().toDouble(),
                this->xElemEdit->text().toInt(), this->yElemEdit->text().toInt(),
                origin, name);
        else
            mesh = MeshLib::MeshGenerator::generateRegularQuadMesh(
                this->xElemEdit->text().toInt(), this->yElemEdit->text().toInt(),
                this->xLengthEdit->text().toDouble(), this->yLengthEdit->text().toDouble(),
                origin, name);
    else if (this->prismButton->isChecked())
        if (this->meshExtentButton->isChecked())
            mesh = MeshLib::MeshGenerator::generateRegularPrismMesh(
                this->xLengthEdit->text().toDouble(), this->yLengthEdit->text().toDouble(),
                this->zLengthEdit->text().toDouble(), this->xElemEdit->text().toInt(),
                this->yElemEdit->text().toInt(), this->zElemEdit->text().toInt(),
                origin, name);
        else
            mesh = MeshLib::MeshGenerator::generateRegularPrismMesh(
                this->xLengthEdit->text().toDouble(), this->yLengthEdit->text().toDouble(),
                this->zLengthEdit->text().toDouble(), this->xElemEdit->text().toInt(),
                this->yElemEdit->text().toInt(), this->zElemEdit->text().toInt(),
                origin, name);
    else if (this->hexButton->isChecked()) {
        if (this->meshExtentButton->isChecked())
        {
            mesh = MeshLib::MeshGenerator::generateRegularHexMesh(
                this->xLengthEdit->text().toDouble(),
                this->yLengthEdit->text().toDouble(),
                this->zLengthEdit->text().toDouble(),
                this->xElemEdit->text().toInt(),
                this->yElemEdit->text().toInt(),
                this->zElemEdit->text().toInt(), origin, name);
        }
        else
        {
            mesh = MeshLib::MeshGenerator::generateRegularHexMesh(
                this->xElemEdit->text().toInt(),
                this->yElemEdit->text().toInt(),
                this->zElemEdit->text().toInt(),
                this->xLengthEdit->text().toDouble(),
                this->yLengthEdit->text().toDouble(),
                this->zLengthEdit->text().toDouble(), origin, name);
        }
    }

    if (mesh == nullptr)
    {
        OGSError::box("Error creating mesh.");
        return;
    }

    auto* const mat_ids = mesh->getProperties().createNewPropertyVector<int>(
        "MaterialIDs", MeshLib::MeshItemType::Cell);
    assert(mat_ids != nullptr);
    mat_ids->resize(mesh->getNumberOfElements());
    emit meshAdded(mesh);
    this->done(QDialog::Accepted);
}

