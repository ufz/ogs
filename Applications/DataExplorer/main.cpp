#include "mainwindow.h"

#include <memory>
#include <QApplication>

#include <logog/include/logog.hpp>
#include "LogogSimpleFormatter.h"
#ifdef VTKFBXCONVERTER_FOUND
#include <fbxsdk.h>
#include "Common.h"
FbxManager* lSdkManager = NULL;
FbxScene* lScene = NULL;
#endif

#include <vtkSmartPointer.h>

#include "BaseLib/BuildInfo.h"
#include "VtkVis/VtkConsoleOutputWindow.h"

int main(int argc, char* argv[])
{
#ifdef VTKFBXCONVERTER_FOUND
    InitializeSdkObjects(lSdkManager, lScene);
#endif

    auto myOutputWindow = vtkSmartPointer<VtkConsoleOutputWindow>::New();
    vtkOutputWindow::SetInstance(myOutputWindow);

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
    std::unique_ptr<MainWindow> w (new MainWindow());
    w->setWindowTitle( w->windowTitle() + " - " +
        QString::fromStdString(BaseLib::BuildInfo::git_describe));
    if (QCoreApplication::arguments().size()>1) {
        w->loadFileOnStartUp(QCoreApplication::arguments().at(1));
    }
    w->show();
    int returncode = a.exec();
    delete formatter;
    delete logogCout;
    LOGOG_SHUTDOWN();
#ifdef VTKFBXCONVERTER_FOUND
    DestroySdkObjects(lSdkManager, true);
#endif

    return returncode;
}
