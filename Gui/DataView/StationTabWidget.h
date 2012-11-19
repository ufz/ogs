/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file StationTabWidget.h
 *
 * Created on 2009-11-03 by Lars Bilke
 */

#ifndef STATIONTABWIDGET_H
#define STATIONTABWIDGET_H

// ** INCLUDES **
#include "ui_StationTabWidgetBase.h"

/**
 * \brief Widget containing StationTreeView-objects.
 */
class StationTabWidget : public QWidget, public Ui_StationTabWidgetBase
{
	Q_OBJECT

public:
	StationTabWidget(QWidget* parent = 0);

private:

private slots:
	void enableSaveButton(bool enable) { this->saveStnPushButton->setEnabled(enable); };
	void enableRemoveButton(bool enable) { this->removeStnPushButton->setEnabled(enable); };
};

#endif // STATIONTABWIDGET_H
