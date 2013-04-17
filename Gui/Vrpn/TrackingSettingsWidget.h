/**
 * \file
 * \author Lars Bilke
 * \date 2010-08-31
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 */

#ifndef TRACKINGSETTINGSWIDGET_H
#define TRACKINGSETTINGSWIDGET_H

#include "ui_TrackingSettingsWidgetBase.h"

class VtkTrackedCamera;
class QSpaceNavigatorClient;
class QVrpnArtTrackingClient;

// TODO Add QSettings stuff
class TrackingSettingsWidget : public QWidget, public Ui_TrackingSettingsWidgetBase
{
	Q_OBJECT

public:
	TrackingSettingsWidget(VtkTrackedCamera* cam, QWidget* parent = 0, Qt::WindowFlags f = 0);
	virtual ~TrackingSettingsWidget();

protected slots:
	void on_deviceNameLineEdit_textChanged() { applyPushButton->setEnabled(true); }
	void on_deviceNameAtLineEdit_textChanged() { applyPushButton->setEnabled(true); }
	void on_applyPushButton_pressed();
	void on_offsetXLineEdit_textChanged(QString text);
	void on_offsetYLineEdit_textChanged(QString text);
	void on_offsetZLineEdit_textChanged(QString text);
	void on_realToVirtualScaleLineEdit_textChanged(QString text);
	void on_aspectRatioLineEdit_textChanged(QString text);
	void on_screenHeightLineEdit_textChanged(QString text);
	void on_nearClipPlaneLineEdit_textChanged(QString text);
	void on_farClipPlaneLineEdit_textChanged(QString text);
private:
	VtkTrackedCamera*  _cam;
};

#endif // TRACKINGSETTINGSWIDGET_H
