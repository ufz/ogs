/**
 * \file MshTabWidget.h
 * 3/11/2009 LB Initial implementation
 * 18/05/2010 KR Re-Implementation
 *
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



};

#endif // MSHTABWIDGET_H
