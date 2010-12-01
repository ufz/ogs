/**
 * \file OsgWidget.h
 * 24/08/2010 LB Initial implementation
 */

#ifndef OSGWIDGET_H
#define OSGWIDGET_H

#include <QGLWidget>
#include <QTimer>

#include <OpenSG/OSGConfig.h>
#include <OpenSG/OSGPassiveWindow.h>
#include <OpenSG/OSGSimpleSceneManager.h>

class OsgWidget : public QGLWidget
{
	Q_OBJECT
	
public:
	OsgWidget(QWidget *parent=0, const QGLWidget *shareWidget=0, Qt::WindowFlags f=0);
	OsgWidget(QGLContext *context, QWidget *parent=0, const QGLWidget *shareWidget=0, Qt::WindowFlags f=0);
	OsgWidget(const QGLFormat &format, QWidget *parent=0, const QGLWidget *shareWidget=0, Qt::WindowFlags f=0);
	virtual ~OsgWidget();
	
	inline OSG::SimpleSceneManager* sceneManager(); 
	inline OSG::PassiveWindowRefPtr window();

protected:
	void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();

private:
	QTimer *_timer;
    OSG::PassiveWindowRefPtr _passiveWin; 
    OSG::SimpleSceneManager* _mgr;
};

inline OSG::SimpleSceneManager* OsgWidget::sceneManager(){
	return _mgr;
}

inline OSG::PassiveWindowRefPtr OsgWidget::window(){
	return _passiveWin;
}

#endif // OSGWIDGET_H
