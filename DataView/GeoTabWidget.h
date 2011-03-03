/**
 * \file GeoTabWidget.h
 * 2011/02/07 KR Initial implementation
 *
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
