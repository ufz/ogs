/**
 * \file StationTabWidget.h
 * 3/11/2009 LB Initial implementation
 *
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

};

#endif // STATIONTABWIDGET_H
