/**
 * \file
 * \author Lars Bilke
 * \date 2010-08-31
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 */

// ** INCLUDES **
#include "QSpaceNavigatorClient.h"
#include "QVrpnArtTrackingClient.h"
#include "TrackingSettingsWidget.h"
#include "VtkTrackedCamera.h"

#include <QCompleter>
#include <QDoubleValidator>
#include <QLineEdit>
#include <QStringList>

TrackingSettingsWidget::TrackingSettingsWidget(VtkTrackedCamera* cam,
                                               QWidget* parent /*= 0*/, Qt::WindowFlags f /*= 0*/)
	: QWidget(parent, f), _cam(cam)
{
	setupUi(this);

	QStringList deviceNameAtSuggestionList;
	deviceNameAtSuggestionList << "141.65.34.36" << "visserv3.intern.ufz.de" << "localhost";
	QCompleter* deviceNameAtCompleter = new QCompleter(deviceNameAtSuggestionList, this);
	deviceNameAtCompleter->setCaseSensitivity(Qt::CaseInsensitive);
	deviceNameAtLineEdit->setCompleter(deviceNameAtCompleter);

	offsetXLineEdit->setValidator(new QDoubleValidator(offsetXLineEdit));
	offsetYLineEdit->setValidator(new QDoubleValidator(offsetYLineEdit));
	offsetZLineEdit->setValidator(new QDoubleValidator(offsetZLineEdit));
	realToVirtualScaleLineEdit->setValidator(new QDoubleValidator(realToVirtualScaleLineEdit));
	aspectRatioLineEdit->setValidator(new QDoubleValidator(aspectRatioLineEdit));
	screenHeightLineEdit->setValidator(new QDoubleValidator(0.1, 10.0, 2, screenHeightLineEdit));
	nearClipPlaneLineEdit->setValidator(new QDoubleValidator(nearClipPlaneLineEdit));
	farClipPlaneLineEdit->setValidator(new QDoubleValidator(farClipPlaneLineEdit));

	QVrpnArtTrackingClient* art = QVrpnArtTrackingClient::Instance();
	QStringList list = art->deviceName().split("@");
	deviceNameLineEdit->setText(list.at(0));
	deviceNameAtLineEdit->setText(list.at(1));
	updateIntervalSpinBox->setValue(art->updateInterval());
	offsetXLineEdit->setText(QString::number(_cam->trackingOffsetX()));
	offsetYLineEdit->setText(QString::number(_cam->trackingOffsetY()));
	offsetZLineEdit->setText(QString::number(_cam->trackingOffsetZ()));
	realToVirtualScaleLineEdit->setText(QString::number(_cam->realToVirtualScale()));
	aspectRatioLineEdit->setText(QString::number(_cam->screenAspectRatio()));
	screenHeightLineEdit->setText(QString::number(_cam->screenHeight()));
}

TrackingSettingsWidget::~TrackingSettingsWidget()
{
}

void TrackingSettingsWidget::on_offsetXLineEdit_textChanged(QString text)
{
	double value = text.toDouble();
	_cam->setTrackingOffsetX(value);
}

void TrackingSettingsWidget::on_offsetYLineEdit_textChanged(QString text)
{
	double value = text.toDouble();
	_cam->setTrackingOffsetY(value);
}
void TrackingSettingsWidget::on_offsetZLineEdit_textChanged(QString text)
{
	double value = text.toDouble();
	_cam->setTrackingOffsetZ(value);
}

void TrackingSettingsWidget::on_realToVirtualScaleLineEdit_textChanged(QString text)
{
	double value = text.toDouble();
	_cam->setRealToVirtualScale(value);
}

void TrackingSettingsWidget::on_aspectRatioLineEdit_textChanged(QString text)
{
	double value = text.toDouble();
	_cam->setScreenAspectRatio(value);
}

void TrackingSettingsWidget::on_screenHeightLineEdit_textChanged(QString text)
{
	double value = text.toDouble();
	_cam->setScreenHeight(value);
}

void TrackingSettingsWidget::on_applyPushButton_pressed()
{
	QVrpnArtTrackingClient* art = QVrpnArtTrackingClient::Instance();
	art->StopTracking();
	QString deviceName = deviceNameLineEdit->text() + "@" + deviceNameAtLineEdit->text();
	art->StartTracking(deviceName.toStdString().c_str(), updateIntervalSpinBox->value());
	applyPushButton->setEnabled(false);
}

void TrackingSettingsWidget::on_nearClipPlaneLineEdit_textChanged( QString text )
{
	double value = text.toDouble();
	double dnear,dfar;
	_cam->GetClippingRange(dnear, dfar);
	_cam->SetClippingRange(value, dfar);
}

void TrackingSettingsWidget::on_farClipPlaneLineEdit_textChanged( QString text )
{
	double value = text.toDouble();
	double dnear, dfar;
	_cam->GetClippingRange(dnear, dfar);
	_cam->SetClippingRange(dnear, value);
}
