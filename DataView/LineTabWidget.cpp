/**
 * \file LineTabWidget.cpp
 * 3/11/2009 LB Initial implementation
 * 
 * Implementation of LineTabWidget
 */

// ** INCLUDES **
#include "LineTabWidget.h"
#include "LinesModel.h"

#include <QMenu>
#include <QContextMenuEvent>


LineTabWidget::LineTabWidget( QWidget* parent /*= 0*/ )
: QWidget(parent)
{
	setupUi(this);
}

void LineTabWidget::contextMenuEvent( QContextMenuEvent* event )
{
	QMenu menu;
	QAction* connectPlyAction = menu.addAction("Connect Polylines...");
	connect(connectPlyAction, SIGNAL(triggered()), static_cast<PolylinesModel*>(this->dataViewWidget->dataView->model()), SLOT(callEditPlyDialog()));
	menu.exec(event->globalPos());
}
