#include "mainwindow.h"

#include <QApplication>
#include <logog/include/logog.hpp>
#include <memory>

#ifdef VTKFBXCONVERTER_FOUND
#include <fbxsdk.h>
#include "ThirdParty/VtkFbxConverter/Common.h"
FbxManager* lSdkManager = nullptr;
FbxScene* lScene = nullptr;
#endif

#include <vtkSmartPointer.h>

#include "BaseLib/BuildInfo.h"
#include "BaseLib/LogogSimpleFormatter.h"
#include "VtkVis/VtkConsoleOutputWindow.h"

int main(int argc, char* argv[])
{
#ifdef VTKFBXCONVERTER_FOUND
    InitializeSdkObjects(lSdkManager, lScene);
#endif

    auto myOutputWindow = vtkSmartPointer<VtkConsoleOutputWindow>::New();
    vtkOutputWindow::SetInstance(myOutputWindow);

    LOGOG_INITIALIZE();
    auto* logogCout = new logog::Cout;
    auto* formatter = new BaseLib::LogogSimpleFormatter;
    logogCout->SetFormatter(*formatter);
    QApplication a(argc, argv);
    QApplication::setApplicationName("OpenGeoSys - Data Explorer");
    QApplication::setApplicationVersion(QString::fromStdString(
        BaseLib::BuildInfo::ogs_version));
    QApplication::setOrganizationName("OpenGeoSys Community");
    QApplication::setOrganizationDomain("opengeosys.org");
    setlocale(LC_NUMERIC,"C");
    QLocale::setDefault(QLocale::German);
    auto w = std::make_unique<MainWindow>();
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
