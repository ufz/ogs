/**
 * \file
 * \author Karsten Rink
 * \date   2010-12-14
 * \brief  Definition of the ModellingTabWidget class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

// ** INCLUDES **
#include "ui_ModellingTabWidgetBase.h"

/**
 * \brief Widget containing the ProcessView.
 */
class ModellingTabWidget : public QWidget, public Ui_ModellingTabWidgetBase
{
    Q_OBJECT

public:
    ModellingTabWidget(QWidget* parent = nullptr);

private slots:

signals:
};
