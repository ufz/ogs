/**
 * \file VisualizationWidget.cpp
 * 3/11/2009 LB Initial implementation
 *
 * Implementation of VisualizationWidget
 */

// ** INCLUDES **
#include "VisualizationWidget.h"
#include "Configure.h"
#include "Point.h"
#include "VtkPickCallback.h"
#include "VtkCustomInteractorStyle.h"
#include "VtkTrackedCamera.h"

#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkCamera.h>
#include <vtkSmartPointer.h>
#include <vtkCellPicker.h>

#include <vtkInteractorStyleSwitch.h>
#include <vtkInteractorStyleRubberBandZoom.h>
#include <vtkMath.h>
#include <vtkCommand.h>
#include <vtkAxesActor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkPNGWriter.h>
#include <vtkWindowToImageFilter.h>
#include <vtkSmartPointer.h>

#include <QSettings>
#include <QFileDialog>
#include <QLineEdit>
#include <QString>
#include <QInputDialog>
#include <QSettings>
#include <QDir>

#ifdef OGS_USE_VRPN
#include "QSpaceNavigatorClient.h"
#include "VtkTrackedCamera.h"
#include <vtkEventQtSlotConnect.h>
#include "QVrpnArtTrackingClient.h"
#include <QTimer>
#endif // OGS_USE_VRPN

VisualizationWidget::VisualizationWidget( QWidget* parent /*= 0*/ )
: QWidget(parent)
{
	this->setupUi(this);

	_interactorStyle = VtkCustomInteractorStyle::New();
	vtkWidget->GetRenderWindow()->GetInteractor()->SetInteractorStyle(_interactorStyle);

	_vtkPickCallback = VtkPickCallback::New();
	vtkSmartPointer<vtkCellPicker> picker = vtkSmartPointer<vtkCellPicker>::New();
	picker->AddObserver(vtkCommand::EndPickEvent, _vtkPickCallback);
	vtkWidget->GetRenderWindow()->GetInteractor()->SetPicker(picker);

	// BUG Render Window conflicts with VREDs render window
#ifndef OGS_VRED_PLUGIN
	vtkRenderWindow* renderWindow = vtkWidget->GetRenderWindow();
	renderWindow->StereoCapableWindowOn();
	renderWindow->SetStereoTypeToCrystalEyes();
	_vtkRender = vtkRenderer::New();
	renderWindow->AddRenderer(_vtkRender);
	_interactorStyle->SetDefaultRenderer(_vtkRender);
#endif // OGS_VRED_PLUGIN                                                                                                               

	QSettings settings("UFZ", "OpenGeoSys-5");

#ifdef OGS_USE_VRPN
	VtkTrackedCamera* cam = new VtkTrackedCamera(this);
	_vtkRender->SetActiveCamera(cam);
	connect( cam, SIGNAL(viewUpdated()), this, SLOT(updateView()) );           

	//QSpaceNavigatorClient* spacenav = QSpaceNavigatorClient::Instance();
	//spacenav->init("spacenav@localhost", 1000 / 15, SpaceNavigatorClient::Z);
	//cam->setFocalPoint(0, 5.0, 0.5);
	//cam->updateView();
	//spacenav->setTranslationFactor(2.0);
	//connect( spacenav, SIGNAL(translated(double, double, double)), cam, SLOT(setTrackingData(double, double, double)) );
	//connect( spacenav, SIGNAL(translated(double, double, double)), cam, SLOT(translate(double, double, double)) );

	QVrpnArtTrackingClient* art = QVrpnArtTrackingClient::Instance(this);
	if (settings.contains("Tracking/artDeviceName"))
	{
		QString deviceName = settings.value("Tracking/artDeviceName").toString();
		QString deviceNameAt = settings.value("Tracking/artDeviceNameAt").toString();
		art->StartTracking(QString(deviceName + "@" + deviceNameAt).toStdString().c_str(),
						   settings.value("Tracking/artUpdateInterval").toInt());
	}
	else
		art->StartTracking("DTrack@141.65.34.36");
	connect( art, SIGNAL(positionUpdated(double, double, double)),
			 cam, SLOT(setTrackingData(double, double, double)) );

	// Connect the vtk event to the qt slot
	_qtConnect = vtkEventQtSlotConnect::New();
	_qtConnect->Connect(vtkWidget->GetRenderWindow()->GetInteractor(),
						vtkCommand::EndInteractionEvent,
						cam,
						SLOT(updatedFromOutside()));

#endif // OGS_USE_VRPN

	_vtkRender->SetBackground(0.0,0.0,0.0);

	// Restore settings
	stereoToolButton->setChecked(settings.value("stereoEnabled").toBool());
	//if (settings.contains("stereoEyeAngle"))
	//	cam->SetEyeAngle(settings.value("stereoEyeAngle").toDouble());
	//else
	//	cam->SetEyeAngle(2.0);
/*
	if (!stereoToolButton->isChecked())
	{
		eyeAngleLabel->setEnabled(false);
		eyeAngleSlider->setEnabled(false);
	}
*/
	//eyeAngleSlider->setValue((int)(_vtkRender->GetActiveCamera()->GetEyeAngle() * 10));
	
	// Create an orientation marker using vtkAxesActor
	vtkSmartPointer<vtkAxesActor> axesActor = vtkSmartPointer<vtkAxesActor>::New();
	vtkOrientationMarkerWidget* markerWidget = vtkOrientationMarkerWidget::New();
	markerWidget->SetOrientationMarker(axesActor);
	//markerWidget->SetViewport(0.0, 0.0, 0.15, 0.3); // size
	markerWidget->SetInteractor(vtkWidget->GetRenderWindow()->GetInteractor());
	markerWidget->EnabledOn();
	markerWidget->InteractiveOff();
}

VisualizationWidget::~VisualizationWidget()
{
	// Write settings
	QSettings settings("UFZ", "OpenGeoSys-5");
	settings.setValue("stereoEnabled", stereoToolButton->isChecked());
	settings.setValue("stereoEyeAngle", _vtkRender->GetActiveCamera()->GetEyeAngle());

	_interactorStyle->deleteLater();
	_vtkPickCallback->deleteLater();
	#ifdef OGS_USE_VRPN
	_qtConnect->Delete();
	#endif // OGS_USE_VRPN
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
	/*
	vtkCamera* camera = _vtkRender->GetActiveCamera();
	double x,y,z;
	camera->GetFocalPoint(x, y, z);
	std::cout << "Focal point: " << x << " " << y << " " << z << std::endl;
	camera->GetPosition(x, y, z);
	std::cout << "Position: " << x << " " << y << " " << z << std::endl;
	camera->GetClippingRange(x, y);
	std::cout << "Clipping range: " << x << " " << y << std::endl << std::endl;
	*/
	
	vtkWidget->GetRenderWindow()->Render();
}

void VisualizationWidget::showAll()
{
	_vtkRender->ResetCamera();
	vtkCamera* cam = _vtkRender->GetActiveCamera();
	double* fp = cam->GetFocalPoint();
	double* p = cam->GetPosition();
	double dist = sqrt(vtkMath::Distance2BetweenPoints(p, fp));
	cam->SetPosition(fp[0], fp[1], fp[2] + dist);
	cam->SetViewUp(0.0, 1.0, 0.0);
	this->updateView();
}

void VisualizationWidget::on_stereoToolButton_toggled( bool checked )
{
	if (checked)
	{
		vtkWidget->GetRenderWindow()->StereoRenderOn();
		//eyeAngleLabel->setEnabled(true);
		//eyeAngleSlider->setEnabled(true);
	}
	else
	{
		vtkWidget->GetRenderWindow()->StereoRenderOff();
		//eyeAngleLabel->setEnabled(false);
		//eyeAngleSlider->setEnabled(false);
	}

	this->updateView();
}
/*
void VisualizationWidget::on_eyeAngleSlider_valueChanged( int value )
{
	Q_UNUSED(value);
	//_vtkRender->GetActiveCamera()->SetEyeAngle(value / 10.0);
	//updateView();
}
*/
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

void VisualizationWidget::on_showAllPushButton_pressed()
{
	this->showAll();
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
	QSettings settings("UFZ", "OpenGeoSys-5");
	QString filename = QFileDialog::getSaveFileName(this, tr("Save screenshot"),
						settings.value("lastScreenshotDir").toString(), "PNG file (*.png)");
	if (filename.count() > 4)
	{
		bool ok;
		int magnification = QInputDialog::getInt(this, tr("Screenshot magnification"),
								tr("Enter a magnification factor for the resulting image."),
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
	writer->SetInput(windowToImageFilter->GetOutput());
	writer->Write();

	this->updateView();
}
