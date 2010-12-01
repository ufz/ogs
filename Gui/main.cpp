#include <QtGui/QApplication>
#include "mainwindow.h"
#include "Configure.h"
#ifdef OGS_USE_OPENSG
#include <OpenSG/OSGBaseFunctions.h>
#endif

int main(int argc, char *argv[])
{
#ifdef OGS_USE_OPENSG
	OSG::osgInit(argc, argv);
#endif

	QApplication a(argc, argv);
	setlocale(LC_NUMERIC,"C");
	MainWindow* w = new MainWindow();
	w->setWindowTitle( w->windowTitle() + " - " + QString(OGS_VERSION) + " - FirstFloor");
	w->show();
	int returncode = a.exec();
	delete w;

#ifdef OGS_USE_OPENSG
	OSG::osgExit();
#endif // OGS_USE_OPENSG

	return returncode;
}
