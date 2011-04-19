/**
 * \file ConditionTabWidget.cpp
 * 2010/12/14 KR Initial implementation
 *
 * Implementation of StationTabWidget
 */

// ** INCLUDES **
#include "ConditionTabWidget.h"
#include "TreeItem.h"
#include "ConditionModel.h"

ConditionTabWidget::ConditionTabWidget( QWidget* parent /*= 0*/ )
: QWidget(parent)
{
	setupUi(this);
}
