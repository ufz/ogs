/**
 * \file SurfaceTabWidget.h
 * 23/04/2010 KR Initial implementation
 *
 */


#ifndef SURFACETABWIDGET_H
#define SURFACETABWIDGET_H

// ** INCLUDES **
#include "ui_SurfaceTabWidgetBase.h"

/**
 * SurfaceTabWidget
 */
class SurfaceTabWidget : public QWidget, public Ui_SurfaceTabWidgetBase
{
	Q_OBJECT

public:
	SurfaceTabWidget(QWidget* parent = 0);


private:

};

#endif // SURFACETABWIDGET_H
