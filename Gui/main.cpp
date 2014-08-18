#include "mainwindow.h"
#include <QtGui/QApplication>
#ifdef VTKOSGCONVERTER_FOUND
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

#include "BaseLib/BuildInfo.h"

int main(int argc, char* argv[])
{
#ifdef VTKOSGCONVERTER_FOUND
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
	QApplication::setApplicationVersion(QString::fromStdString(
		BaseLib::BuildInfo::ogs_version));
	QApplication::setOrganizationName("OpenGeoSys Community");
	QApplication::setOrganizationDomain("opengeosys.org");
	setlocale(LC_NUMERIC,"C");
	QLocale::setDefault(QLocale::German);
	MainWindow* w = new MainWindow();
	w->setWindowTitle( w->windowTitle() + " - " +
		QString::fromStdString(BaseLib::BuildInfo::ogs_version_and_persons) +
		" - FirstFloor");
	w->show();
	int returncode = a.exec();
	delete w;
	delete formatter;
	delete logogCout;
	LOGOG_SHUTDOWN();
#ifdef VTKFBXCONVERTER_FOUND
	DestroySdkObjects(lSdkManager, true);
#endif
#ifdef VTKOSGCONVERTER_FOUND
	OSG::osgExit();
#endif

	return returncode;
}
