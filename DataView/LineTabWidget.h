/**
 * \file LineTabWidget.h
 * 3/11/2009 LB Initial implementation
 *
 */


#ifndef LINETABWIDGET_H
#define LINETABWIDGET_H

// ** INCLUDES **
#include "ui_LineTabWidgetBase.h"

/**
 * LineTabWidget
 */
class LineTabWidget : public QWidget, public Ui_LineTabWidgetBase
{
	Q_OBJECT

public:
	LineTabWidget(QWidget* parent = 0);


private:

};

#endif // LINETABWIDGET_H
