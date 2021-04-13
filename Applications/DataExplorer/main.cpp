#include "mainwindow.h"

#include <QApplication>
#include <QSurfaceFormat>
#include <QVTKOpenGLWidget.h>
#include <memory>
#include <vtkSmartPointer.h>

#include "InfoLib/GitInfo.h"
#include "VtkVis/VtkConsoleOutputWindow.h"

// TODO: Replace this on VTK 9 upgrade, see
// https://discourse.vtk.org/t/vtk-use-file/3645/2
#if VTK_VIA_CPM
#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkInteractionStyle)
VTK_MODULE_INIT(vtkRenderingFreeType)
VTK_MODULE_INIT(vtkRenderingOpenGL2)
#endif

int main(int argc, char* argv[])
{
    // needed to ensure appropriate OpenGL context is created for VTK rendering.
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLWidget::defaultFormat());

    auto myOutputWindow = vtkSmartPointer<VtkConsoleOutputWindow>::New();
    vtkOutputWindow::SetInstance(myOutputWindow);

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

    return returncode;
}
