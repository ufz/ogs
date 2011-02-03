/**
 * \file PntTabWidget.h
 * 3/11/2009 LB Initial implementation
 *
 */


#ifndef PNTTABWIDGET_H
#define PNTTABWIDGET_H

// ** INCLUDES **
#include "ui_PntTabWidgetBase.h"

/**
 * \brief Widget containing the view for the PointsModel.
 */
class PntTabWidget : public QWidget, public Ui_PntTabWidgetBase
{
	Q_OBJECT

public:
	PntTabWidget(QWidget* parent = 0);

};

#endif // PNTTABWIDGET_H
