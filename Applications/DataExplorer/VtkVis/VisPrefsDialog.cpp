/**
 * \file
 * \author Karsten Rink
 * \date   2010-06-14
 * \brief  Implementation of the VisPrefsDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "VisPrefsDialog.h"
#include <QDoubleValidator>
#include <QLineEdit>
#include <QSettings>
#include <QVariant>

#include "VisualizationWidget.h"
#include "VtkVisPipeline.h"

/// Constructor
VisPrefsDialog::VisPrefsDialog(VtkVisPipeline &pipeline,
                               VisualizationWidget &widget,
                               QDialog* parent) :
    QDialog(parent), vtkVisPipeline_(pipeline), visWidget_(widget),
    above_(0,0,2000000), below_(0,0,-2000000)
{
    setupUi(this);
    if (vtkVisPipeline_.getLight(above_))
    {
        lightAboveBox->toggle();
    }
    if (vtkVisPipeline_.getLight(below_))
    {
        lightBelowBox->toggle();
    }

    bgColorButton->setColor(vtkVisPipeline_.getBGColor());

    QValidator* validator = new QDoubleValidator(0, 100000, 2, this);
    superelevationLineEdit->setValidator(validator);
    QSettings settings;
    superelevationLineEdit->setText(settings.value("globalSuperelevation", 1.0).toString());
    cullBackfacesCheckBox->setCheckState(Qt::CheckState(settings.value("globalCullBackfaces", 0).toInt()));
    loadShowAllCheckBox->setCheckState(Qt::CheckState(settings.value("resetViewOnLoad", 2).toInt()));
}

void VisPrefsDialog::on_bgColorButton_colorPicked( QColor color )
{
    QColor bgColor(color.red(), color.green(), color.blue());
    vtkVisPipeline_.setBGColor(bgColor);
}

void VisPrefsDialog::on_lightAboveBox_clicked()
{
    if (lightAboveBox->isChecked())
    {
        vtkVisPipeline_.addLight(above_);
    }
    else
    {
        vtkVisPipeline_.removeLight(above_);
    }
}

void VisPrefsDialog::on_lightBelowBox_clicked()
{
    if (lightBelowBox->isChecked())
    {
        vtkVisPipeline_.addLight(below_);
    }
    else
    {
        vtkVisPipeline_.removeLight(below_);
    }
}

void VisPrefsDialog::on_superelevationPushButton_pressed()
{
    double factor = superelevationLineEdit->text().toDouble();
    vtkVisPipeline_.setGlobalSuperelevation(factor);

    QSettings settings;
    settings.setValue("globalSuperelevation", factor);
}

void VisPrefsDialog::on_loadShowAllCheckBox_stateChanged(int state)
{
    visWidget_.setShowAllOnLoad(static_cast<bool>(state));
    vtkVisPipeline_.resetCameraOnAddOrRemove(static_cast<bool>(state));

    QSettings settings;
    settings.setValue("resetViewOnLoad", state);
}

void VisPrefsDialog::on_cullBackfacesCheckBox_stateChanged(int state)
{
    vtkVisPipeline_.setGlobalBackfaceCulling(static_cast<bool>(state));

    QSettings settings;
    settings.setValue("globalCullBackfaces", state);
}

