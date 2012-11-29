/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file MshTabWidget.h
 *
 * Created on 2009-11-03 by Lars Bilke
 */

#ifndef MSHTABWIDGET_H
#define MSHTABWIDGET_H

// ** INCLUDES **
#include "ui_MshTabWidgetBase.h"

/**
 * \brief Widget for data views of meshes.
 */
class MshTabWidget : public QWidget, public Ui_MshTabWidgetBase
{
	Q_OBJECT

public:
	MshTabWidget(QWidget* parent = 0);

private slots:
	void enableSaveButton(bool enable) { this->saveMeshPushButton->setEnabled(enable); };
	void enableRemoveButton(bool enable) { this->removeMeshPushButton->setEnabled(enable); };
};

#endif // MSHTABWIDGET_H
