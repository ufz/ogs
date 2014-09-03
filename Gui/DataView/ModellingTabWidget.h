/**
 * \file
 * \author Karsten Rink
 * \date   2010-12-14
 * \brief  Definition of the ModellingTabWidget class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MODELLINGTABWIDGET_H
#define MODELLINGTABWIDGET_H

// ** INCLUDES **
#include "ui_ModellingTabWidgetBase.h"

/**
 * \brief Widget containing the ProcessView.
 */
class ModellingTabWidget : public QWidget, public Ui_ModellingTabWidgetBase
{
	Q_OBJECT

public:
	ModellingTabWidget(QWidget* parent = 0);

private slots:
	void on_addProcessButton_pressed();
	void on_deleteAllButton_pressed();

signals:
	void requestNewProcess();
};

#endif // MODELLINGTABWIDGET_H
