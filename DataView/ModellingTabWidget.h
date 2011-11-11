/**
 * \file ModellingTabWidget.h
 * 2010/12/14 KR Initial implementation
 *
 */

#ifndef MODELLINGTABWIDGET_H
#define MODELLINGTABWIDGET_H

// ** INCLUDES **
#include "ui_ConditionTabWidgetBase.h"

/**
 * \brief Widget containing the ProcessView.
 */
class ModellingTabWidget : public QWidget, public Ui_ConditionTabWidgetBase
{
	Q_OBJECT

public:
	ModellingTabWidget(QWidget* parent = 0);
};

#endif // MODELLINGTABWIDGET_H
