/**
 * \file
 * \author Karsten Rink
 * \date   2010-12-14
 * \brief  Implementation of the ModellingTabWidget class.
 *
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ModellingTabWidget.h"

ModellingTabWidget::ModellingTabWidget(QWidget* parent /*= 0*/)
    : QWidget(parent)
{
    setupUi(this);
}
