/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file GeoTabWidget.h
 *
 * Created on 2011-02-07 by Karsten Rink
 */

#ifndef GEOTABWIDGET_H
#define GEOTABWIDGET_H

// ** INCLUDES **
#include "ui_GeoTabWidgetBase.h"

/**
 * \brief Widget containing GeoTreeView-objects.
 */
class GeoTabWidget : public QWidget, public Ui_GeoTabWidgetBase
{
	Q_OBJECT

public:
	GeoTabWidget(QWidget* parent = 0);

private:
};

#endif // GEOTABWIDGET_H
