#include "Configure.h"
#include "mainwindow.h"
#include <QtGui/QApplication>
#ifdef OGS_USE_OPENSG
#include <OpenSG/OSGBaseFunctions.h>
#endif
#include "logog/include/logog.hpp"
#include "LogogSimpleFormatter.h"
#ifdef VTKFBXCONVERTER_FOUND
#include <fbxsdk.h>
#include "Common.h"
FbxManager* lSdkManager = NULL;
FbxScene* lScene = NULL;
#endif

int main(int argc, char* argv[])
{
#ifdef OGS_USE_OPENSG
	OSG::osgInit(argc, argv);
#endif
#ifdef VTKFBXCONVERTER_FOUND
	InitializeSdkObjects(lSdkManager, lScene);
#endif
	LOGOG_INITIALIZE();
	logog::Cout* logogCout = new logog::Cout;
	BaseLib::LogogSimpleFormatter* formatter = new BaseLib::LogogSimpleFormatter;
	logogCout->SetFormatter(*formatter);
	QApplication a(argc, argv);
	QApplication::setApplicationName("OpenGeoSys - Data Explorer");
	QApplication::setApplicationVersion(QString(OGS_VERSION));
	QApplication::setOrganizationName("OpenGeoSys Community");
	QApplication::setOrganizationDomain("opengeosys.org");
	setlocale(LC_NUMERIC,"C");
	MainWindow* w = new MainWindow();
	w->setWindowTitle( w->windowTitle() + " - " + QString(OGS_VERSION_AND_PERSONS) + " - FirstFloor");
	w->show();
	int returncode = a.exec();
	delete w;
	delete formatter;
	delete logogCout;
	LOGOG_SHUTDOWN();
#ifdef VTKFBXCONVERTER_FOUND
	DestroySdkObjects(lSdkManager);
#endif
#ifdef OGS_USE_OPENSG
	OSG::osgExit();
#endif // OGS_USE_OPENSG

	return returncode;
}
