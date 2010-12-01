/**
 * \file VisPrefsDialog.cpp
 * 14/06/2010  KR Initial implementation
 */

#include <QSettings>
#include "VisPrefsDialog.h"

#include "VtkVisPipeline.h"

/// Constructor
VisPrefsDialog::VisPrefsDialog(VtkVisPipeline* pipeline, QDialog* parent) : 
	QDialog(parent), _vtkVisPipeline(pipeline), _above(0,0,2000000), _below(0,0,-2000000)
{
	setupUi(this);
	if (_vtkVisPipeline->getLight(_above))
		lightAboveBox->toggle();
	if (_vtkVisPipeline->getLight(_below))
		lightBelowBox->toggle();

	QColor color = _vtkVisPipeline->getBGColor();
	bgColorButton->setColor(_vtkVisPipeline->getBGColor());
}

VisPrefsDialog::~VisPrefsDialog()
{
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

void VisPrefsDialog::accept()
{
	this->done(QDialog::Accepted);
}


void VisPrefsDialog::reject()
{
	this->done(QDialog::Rejected);
}
