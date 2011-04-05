/**
 * \file ConditionTabWidget.h
 * 2010/12/14 KR Initial implementation
 *
 */


#ifndef CONDITIONTABWIDGET_H
#define CONDITIONTABWIDGET_H

// ** INCLUDES **
#include "ui_ConditionTabWidgetBase.h"

/**
 * \brief Widget containing the ConditionView.
 */
class ConditionTabWidget : public QWidget, public Ui_ConditionTabWidgetBase
{
	Q_OBJECT

public:
	ConditionTabWidget(QWidget* parent = 0);

};

#endif // CONDITIONTABWIDGET_H
