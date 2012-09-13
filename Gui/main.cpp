#include "Configure.h"
#include "mainwindow.h"
#include <QtGui/QApplication>
#ifdef OGS_USE_OPENSG
#include <OpenSG/OSGBaseFunctions.h>
#endif
#include "logog.hpp"
#include "LogogSimpleFormatter.h"

int main(int argc, char* argv[])
{
#ifdef OGS_USE_OPENSG
	OSG::osgInit(argc, argv);
#endif
	LOGOG_INITIALIZE();
	logog::Cout* logogCout = new logog::Cout;
	BaseLib::LogogSimpleFormatter formatter;
	logogCout->SetFormatter(formatter);
	QApplication a(argc, argv);
	setlocale(LC_NUMERIC,"C");
	MainWindow* w = new MainWindow();
	w->setWindowTitle( w->windowTitle() + " - " + QString(OGS_VERSION) + " - FirstFloor");
	w->show();
	int returncode = a.exec();
	delete w;

	delete logogCout;
	LOGOG_SHUTDOWN();
#ifdef OGS_USE_OPENSG
	OSG::osgExit();
#endif // OGS_USE_OPENSG

	return returncode;
}
