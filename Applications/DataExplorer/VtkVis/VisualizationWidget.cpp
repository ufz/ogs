/**
 * \file
 * \author Lars Bilke
 * \date   2009-11-03
 * \brief  Implementation of the VisualizationWidget class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include <cmath>

#include "Point.h"
#include "VisualizationWidget.h"
#include "VtkCustomInteractorStyle.h"
#include "VtkPickCallback.h"

#include <vtkCamera.h>
#include <vtkCellPicker.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>

#include <vtkAxesActor.h>
#include <vtkCommand.h>
#include <vtkInteractorStyleRubberBandZoom.h>
#include <vtkInteractorStyleSwitch.h>
#include <vtkMath.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkPNGWriter.h>
#include <vtkSmartPointer.h>
#include <vtkWindowToImageFilter.h>

#include <QCursor>
#include <QDir>
#include <QFileDialog>
#include <QInputDialog>
#include <QLineEdit>
#include <QSettings>
#include <QString>

VisualizationWidget::VisualizationWidget( QWidget* parent /*= 0*/ )
: QWidget(parent), _vtkRender(nullptr), _markerWidget(nullptr),
  _interactorStyle(nullptr), _vtkPickCallback(nullptr)
{
    this->setupUi(this);

    _interactorStyle = VtkCustomInteractorStyle::New();
    vtkWidget->GetRenderWindow()->GetInteractor()->SetInteractorStyle(_interactorStyle);

    _vtkPickCallback = VtkPickCallback::New();
    vtkSmartPointer<vtkCellPicker> picker = vtkSmartPointer<vtkCellPicker>::New();
    picker->AddObserver(vtkCommand::EndPickEvent, _vtkPickCallback);
    vtkWidget->GetRenderWindow()->GetInteractor()->SetPicker(picker);

    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkWidget->GetRenderWindow();
    renderWindow->StereoCapableWindowOn();
    renderWindow->SetStereoTypeToCrystalEyes();
    _vtkRender = vtkRenderer::New();
    renderWindow->AddRenderer(_vtkRender);
    _interactorStyle->SetDefaultRenderer(_vtkRender);

    QSettings settings;

    _vtkRender->SetBackground(0.0,0.0,0.0);

    // Create an orientation marker using vtkAxesActor
    vtkSmartPointer<vtkAxesActor> axesActor = vtkSmartPointer<vtkAxesActor>::New();
    _markerWidget = vtkOrientationMarkerWidget::New();
    _markerWidget->SetOrientationMarker(axesActor);
    _markerWidget->SetInteractor(vtkWidget->GetRenderWindow()->GetInteractor());
    _markerWidget->EnabledOn();
    _markerWidget->InteractiveOff();

    _isShowAllOnLoad = settings.value("resetViewOnLoad", true).toBool();

    // Set alternate cursor shapes
    connect(_interactorStyle, SIGNAL(cursorChanged(Qt::CursorShape)),
            this, SLOT(setCursorShape(Qt::CursorShape)));
}

VisualizationWidget::~VisualizationWidget()
{
    _vtkPickCallback->Delete();
    _interactorStyle->Delete();
    _markerWidget->Delete();
    _vtkRender->Delete();
}

VtkCustomInteractorStyle* VisualizationWidget::interactorStyle() const
{
    return _interactorStyle;
}

VtkPickCallback* VisualizationWidget::vtkPickCallback() const
{
    return _vtkPickCallback;
}

void VisualizationWidget::updateView()
{
    if(vtkWidget->GetRenderWindow()->IsDrawable())
        vtkWidget->GetRenderWindow()->Render();
}

void VisualizationWidget::showAll(int x, int y, int z)
{
    _vtkRender->ResetCamera();
    vtkCamera* cam = _vtkRender->GetActiveCamera();
    double* fp = cam->GetFocalPoint();
    double* p = cam->GetPosition();
    double dist = sqrt(vtkMath::Distance2BetweenPoints(p, fp));
    cam->SetPosition(fp[0]+(x*dist), fp[1]+(y*dist), fp[2]+(z*dist));

    if (x!=0 || y!=0)
        cam->SetViewUp(0.0, 0.0, 1.0);
    else
        cam->SetViewUp(0.0, 1.0, 0.0);
    this->updateView();
}

void VisualizationWidget::updateViewOnLoad()
{
    if (_isShowAllOnLoad)
        this->showAll(0, 0, 1);
    else
        updateView();
}

void VisualizationWidget::on_zoomToolButton_toggled( bool checked )
{
    if (checked)
    {
        vtkSmartPointer<vtkInteractorStyleRubberBandZoom> interactorStyle =
                vtkSmartPointer<vtkInteractorStyleRubberBandZoom>::New();
        vtkWidget->GetRenderWindow()->GetInteractor()->SetInteractorStyle(interactorStyle);
        QCursor cursor;
        cursor.setShape(Qt::CrossCursor);
        vtkWidget->setCursor(cursor);
    }
    else
    {
        vtkWidget->GetRenderWindow()->GetInteractor()->SetInteractorStyle(_interactorStyle);
        QCursor cursor;
        cursor.setShape(Qt::ArrowCursor);
        vtkWidget->setCursor(cursor);
    }
}

void VisualizationWidget::on_highlightToolButton_toggled(bool checked)
{
    _interactorStyle->setHighlightActor(checked);
}

void VisualizationWidget::on_orthogonalProjectionToolButton_toggled( bool checked )
{
    _vtkRender->GetActiveCamera()->SetParallelProjection(checked);
    this->updateView();
}

void VisualizationWidget::on_screenshotPushButton_pressed()
{
    QSettings settings;
    QString filename = QFileDialog::getSaveFileName(this, tr("Save screenshot"),
                                                    settings.value(
                                                            "lastScreenshotDir").toString(),
                                                    "PNG file (*.png)");
    if (filename.count() > 4)
    {
        bool ok;
        int magnification = QInputDialog::getInt(this, tr("Screenshot magnification"),
                                                 tr(
                                                         "Enter a magnification factor for the resulting image."),
                                                 2, 1, 10, 1, &ok);
        if (ok)
        {
            QDir dir(filename);
            settings.setValue("lastScreenshotDir", dir.absolutePath());
            this->screenshot(filename, magnification);
        }
    }
}

void VisualizationWidget::screenshot(QString filename, int magnification)
{
    vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter =
        vtkSmartPointer<vtkWindowToImageFilter>::New();
    windowToImageFilter->SetInput(vtkWidget->GetRenderWindow());
    // Set the resolution of the output image
    // magnification times the current resolution of vtk render window
    windowToImageFilter->SetMagnification(magnification);
    // Also record the alpha (transparency) channel
    windowToImageFilter->SetInputBufferTypeToRGBA();
    windowToImageFilter->Update();

    vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
    writer->SetFileName(filename.toStdString().c_str());
    writer->SetInputData(windowToImageFilter->GetOutput());
    writer->Write();

    this->updateView();
}

void VisualizationWidget::setCursorShape(Qt::CursorShape shape)
{
    this->setCursor(QCursor(shape));
}
