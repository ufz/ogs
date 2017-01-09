/**
 * \file
 * \author Karsten Rink
 * \date   2012-01-04
 * \brief  Implementation of the CondFromRasterDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CondFromRasterDialog.h"
#include "Mesh.h"

#include <QFileDialog>
#include <QSettings>

#include "DirectConditionGenerator.h"
#include "OGSError.h"
#include "StrictDoubleValidator.h"

CondFromRasterDialog::CondFromRasterDialog(const std::vector<MeshLib::Mesh*> &msh_vec, QDialog* parent)
    : QDialog(parent), _msh_vec(msh_vec)
{
    setupUi(this);

    this->scalingEdit->setEnabled(false);
    _scale_validator = new StrictDoubleValidator(-1e+10, 1e+20, 5);
    this->scalingEdit->setText("1.0");
    this->scalingEdit->setValidator (_scale_validator);

    for (std::vector<MeshLib::Mesh*>::const_iterator it = _msh_vec.begin(); it != _msh_vec.end(); ++it)
        this->meshBox->addItem(QString::fromStdString((*it)->getName()));

    this->directButton->setChecked(true);
}

CondFromRasterDialog::~CondFromRasterDialog()
{
    delete _scale_validator;
}

void CondFromRasterDialog::on_selectButton_pressed()
{
    QSettings settings;
#ifdef GEOTIFF_FOUND
    QString geotiffExtension(" *.tif");
#else
    QString geotiffExtension("");
#endif
    QString fileName = QFileDialog::getOpenFileName(this, "Select raster file",
                                                    settings.value("lastOpenedRasterFileDirectory").toString(),
                                                    QString("Raster files (*.asc *.grd);;").arg(geotiffExtension));

    if (!fileName.isEmpty())
    {
        this->rasterEdit->setText(fileName);

        QFileInfo fi(fileName);
        settings.setValue("lastOpenedRasterFileDirectory", fi.absolutePath());
    }
}


void CondFromRasterDialog::accept()
{
    std::string mesh_name (this->meshBox->currentText().toStdString());
    std::string raster_name (this->rasterEdit->text().toStdString());
    double scaling_factor = this->scalingEdit->text().toDouble();
    std::vector< std::pair<std::size_t,double> > direct_values;

    if (mesh_name.empty())
    {
        OGSError::box("No mesh selected.");
        return;
    }
    if (raster_name.empty())
    {
        OGSError::box("No raster selected.");
        return;
    }

    MeshLib::Mesh* mesh (NULL);
    for (std::size_t i=0; i<_msh_vec.size(); i++)
        if (_msh_vec[i]->getName().compare(mesh_name) == 0)
        {
            mesh = _msh_vec[i];
            break;
        }


    if (this->directButton->isChecked())
    {
        DirectConditionGenerator dcg;
        direct_values = dcg.directToSurfaceNodes(*mesh, raster_name);
        //dcg.writeToFile(direct_node_name);
    }
    else
    {
        if (scaling_factor <= 0)
        {
            OGSError::box("No valid scaling factor given.");
            return;
        }
        MeshLib::Mesh* new_mesh = const_cast<MeshLib::Mesh*>(mesh);
        DirectConditionGenerator dcg;
        direct_values = dcg.directWithSurfaceIntegration(*new_mesh, raster_name, scaling_factor);

        //dcg.writeToFile(direct_node_name);
    }
    //emit directNodesWritten(direct_node_name);
    emit transmitDisValues(direct_values);
    this->done(QDialog::Accepted);
}

void CondFromRasterDialog::reject()
{
    this->done(QDialog::Rejected);
}

void CondFromRasterDialog::on_integrateButton_toggled(bool isSelected)
{
    this->scalingEdit->setEnabled(isSelected);
}
