/**
 * \file
 * \author Karsten Rink
 * \date   2010-06-14
 * \brief  Implementation of the VisPrefsDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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
    QDialog(parent), _vtkVisPipeline(pipeline), _visWidget(widget),
    _above(0,0,2000000), _below(0,0,-2000000)
{
    setupUi(this);
    if (_vtkVisPipeline.getLight(_above))
        lightAboveBox->toggle();
    if (_vtkVisPipeline.getLight(_below))
        lightBelowBox->toggle();

    bgColorButton->setColor(_vtkVisPipeline.getBGColor());

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
    _vtkVisPipeline.setBGColor(bgColor);
}

void VisPrefsDialog::on_lightAboveBox_clicked()
{
    if (lightAboveBox->isChecked())
        _vtkVisPipeline.addLight(_above);
    else
        _vtkVisPipeline.removeLight(_above);
}

void VisPrefsDialog::on_lightBelowBox_clicked()
{
    if (lightBelowBox->isChecked())
        _vtkVisPipeline.addLight(_below);
    else
        _vtkVisPipeline.removeLight(_below);
}

void VisPrefsDialog::on_superelevationPushButton_pressed()
{
    double factor = superelevationLineEdit->text().toDouble();
    _vtkVisPipeline.setGlobalSuperelevation(factor);

    QSettings settings;
    settings.setValue("globalSuperelevation", factor);
}

void VisPrefsDialog::on_loadShowAllCheckBox_stateChanged(int state)
{
    _visWidget.setShowAllOnLoad((bool)state);
    _vtkVisPipeline.resetCameraOnAddOrRemove((bool)state);

    QSettings settings;
    settings.setValue("resetViewOnLoad", state);
}

void VisPrefsDialog::on_cullBackfacesCheckBox_stateChanged(int state)
{
    _vtkVisPipeline.setGlobalBackfaceCulling((bool)state);

    QSettings settings;
    settings.setValue("globalCullBackfaces", state);
}

