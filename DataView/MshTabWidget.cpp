/**
 * \file MshTabWidget.cpp
 * 3/11/2009 LB Initial implementation
 * 18/05/2010 KR Re-Implementation
 *
 * Implementation of MshTabWidget
 */

// ** INCLUDES **
#include "MshTabWidget.h"

MshTabWidget::MshTabWidget( QWidget* parent /*= 0*/ )
	: QWidget(parent)
{
	setupUi(this);

	connect(this->addMeshPushButton, SIGNAL(clicked()), this->treeView, SLOT(addMeshAction()));
	connect(this->clearAllPushButton, SIGNAL(clicked()), this->treeView, SLOT(removeAllMeshes()));
}
