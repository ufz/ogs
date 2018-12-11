#include "mainwindow.h"

#include <QApplication>
#include <QSurfaceFormat>
#include <QVTKOpenGLWidget.h>
#include <memory>
#include "BaseLib/Logging.h"

#ifdef VTKFBXCONVERTER_FOUND
#include <fbxsdk.h>
#include "ThirdParty/VtkFbxConverter/Common.h"
FbxManager* lSdkManager = nullptr;
FbxScene* lScene = nullptr;
#endif

#include <vtkSmartPointer.h>

#include "InfoLib/GitInfo.h"
#include "BaseLib/LogogSimpleFormatter.h"
#include "VtkVis/VtkConsoleOutputWindow.h"

int main(int argc, char* argv[])
{
#ifdef VTKFBXCONVERTER_FOUND
    InitializeSdkObjects(lSdkManager, lScene);
#endif

    // needed to ensure appropriate OpenGL context is created for VTK rendering.
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLWidget::defaultFormat());

    auto myOutputWindow = vtkSmartPointer<VtkConsoleOutputWindow>::New();
    vtkOutputWindow::SetInstance(myOutputWindow);

    LOGOG_INITIALIZE();
    auto* logogCout = new logog::Cout;
    auto* formatter = new BaseLib::LogogSimpleFormatter;
    logogCout->SetFormatter(*formatter);
    QApplication a(argc, argv);
    QApplication::setApplicationName("OpenGeoSys - Data Explorer");
    QApplication::setApplicationVersion(QString::fromStdString(
        GitInfoLib::GitInfo::ogs_version));
    QApplication::setOrganizationName("OpenGeoSys Community");
    QApplication::setOrganizationDomain("opengeosys.org");
    setlocale(LC_NUMERIC,"C");
    QLocale::setDefault(QLocale::German);
    auto w = std::make_unique<MainWindow>();
    w->setWindowTitle( w->windowTitle() + " - " +
        QString::fromStdString(GitInfoLib::GitInfo::ogs_version));
    if (QCoreApplication::arguments().size()>1) {
        w->loadFileOnStartUp(QCoreApplication::arguments().at(1));
    }
    w->show();
    int returncode = QApplication::exec();
    delete formatter;
    delete logogCout;
    LOGOG_SHUTDOWN();
#ifdef VTKFBXCONVERTER_FOUND
    DestroySdkObjects(lSdkManager, true);
#endif

    return returncode;
}
