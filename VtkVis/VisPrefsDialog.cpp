/**
 * \file VisPrefsDialog.cpp
 * 14/06/2010  KR Initial implementation
 */

#include <QSettings>
#include <QDoubleValidator>
#include <QLineEdit>
#include <QVariant>
#include "VisPrefsDialog.h"

#include "VtkVisPipeline.h"
#include "VisualizationWidget.h"

/// Constructor
VisPrefsDialog::VisPrefsDialog(VtkVisPipeline* pipeline, VisualizationWidget* widget, QDialog* parent) :
	QDialog(parent), _vtkVisPipeline(pipeline), _visWidget(widget), _above(0,0,2000000), _below(0,0,-2000000)
{
	setupUi(this);
	if (_vtkVisPipeline->getLight(_above))
		lightAboveBox->toggle();
	if (_vtkVisPipeline->getLight(_below))
		lightBelowBox->toggle();

	QColor color = _vtkVisPipeline->getBGColor();
	bgColorButton->setColor(_vtkVisPipeline->getBGColor());

	QValidator* validator = new QDoubleValidator(0, 100000, 2, this);
	superelevationLineEdit->setValidator(validator);
	QSettings settings("UFZ, OpenGeoSys-5");
	superelevationLineEdit->setText(settings.value("globalSuperelevation", 1.0).toString());
}

void VisPrefsDialog::on_bgColorButton_colorPicked( QColor color )
{
	QColor bgColor(color.red(), color.green(), color.blue());
	_vtkVisPipeline->setBGColor(bgColor);
}

void VisPrefsDialog::on_lightAboveBox_clicked()
{
	if (lightAboveBox->isChecked())
		_vtkVisPipeline->addLight(_above);
	else
		_vtkVisPipeline->removeLight(_above);
}

void VisPrefsDialog::on_lightBelowBox_clicked()
{
	if (lightBelowBox->isChecked())
		_vtkVisPipeline->addLight(_below);
	else
		_vtkVisPipeline->removeLight(_below);
}

void VisPrefsDialog::on_superelevationPushButton_pressed()
{
	double factor = superelevationLineEdit->text().toDouble();
	_vtkVisPipeline->setGlobalSuperelevation(factor);

	QSettings settings("UFZ, OpenGeoSys-5");
	settings.setValue("globalSuperelevation", factor);
}

void VisPrefsDialog::on_loadShowAllCheckBox_stateChanged(int state)
{
	_visWidget->setShowAllOnLoad((bool)state);
	_vtkVisPipeline->resetCameraOnAddOrRemove((bool)state);
}
