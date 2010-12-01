/**
 * \file OsgWidget.cpp
 * 24/08/2010 LB Initial implementation
 * 
 * Implementation of OsgWidget class
 */

// ** INCLUDES **
#include "OsgWidget.h"

OsgWidget::OsgWidget (QWidget * parent, const QGLWidget *shareWidget, Qt::WindowFlags f)
:QGLWidget(parent, shareWidget, f)
{
	_passiveWin = OSG::PassiveWindow::create();					 
	_mgr = new OSG::SimpleSceneManager();
	_mgr->setWindow(_passiveWin);
}

OsgWidget::OsgWidget (QGLContext *context, QWidget *parent, const QGLWidget *shareWidget, Qt::WindowFlags f)
:QGLWidget(context, parent, shareWidget, f)
{
   	 _passiveWin = OSG::PassiveWindow::create();					 
	_mgr = new OSG::SimpleSceneManager();
	_mgr->setWindow(_passiveWin); 
}

OsgWidget::OsgWidget (const QGLFormat &format, QWidget *parent, const QGLWidget *shareWidget, Qt::WindowFlags f)
:QGLWidget(format, parent, shareWidget, f)
{
   	 _passiveWin = OSG::PassiveWindow::create();					 
	_mgr = new OSG::SimpleSceneManager();
	_mgr->setWindow(_passiveWin);
}

OsgWidget::~OsgWidget()
{
	delete _mgr;
}

void OsgWidget::initializeGL()
{
   	 _passiveWin->init();
   	 _timer = new QTimer(this);
   	 connect(_timer, SIGNAL(timeout()), this, SLOT(updateGL()));
   	 _timer->start(2000);			 
}

void OsgWidget::resizeGL(int w, int h)
{
   	 _passiveWin->resize(w,h);
}

void OsgWidget::paintGL()
{
    _mgr->redraw();
}