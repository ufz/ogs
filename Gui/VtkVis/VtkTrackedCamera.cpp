/**
 * \file VtkTrackedCamera.cpp
 * 25/08/2010 LB Initial implementation
 *
 * Implementation of VtkTrackedCamera class
 */

// ** INCLUDES **
#include "VtkTrackedCamera.h"

#include <vtkMath.h>
#include <vtkPerspectiveTransform.h>

#include <QSettings>

VtkTrackedCamera::VtkTrackedCamera(QObject* parent)
	: QObject(parent), vtkOpenGLCamera()
{
	QSettings settings("UFZ", "OpenGeoSys-5");
	settings.beginGroup("Tracking");
	// Read from settings
	if (settings.contains("artOffsetX"))
	{
		_trackedPositionOffset[0] = settings.value("artOffsetX").toDouble();
		_trackedPositionOffset[1] = settings.value("artOffsetY").toDouble();
		_trackedPositionOffset[2] = settings.value("artOffsetZ").toDouble();
		_realToVirtualScale = settings.value("artRealToVirtualScale").toDouble();
		_screenAspectRatio = settings.value("artAspectRatio").toDouble();
		_screenHeight = settings.value("artScreenHeight").toDouble();
	}
	// Default values
	else
	{
		_trackedPositionOffset[0] = 0;
		_trackedPositionOffset[1] = 0;
		_trackedPositionOffset[2] = 0;
		_realToVirtualScale = 1.0;
		_screenAspectRatio = 1.6;
		_screenHeight = 0.4;
	}
	settings.endGroup();

	// if z up
	_trackedPosition[0] = 0;
	_trackedPosition[1] = 2;
	_trackedPosition[2] = 0;
	_focalPoint[0] = 0;
	_focalPoint[1] = 0;
	_focalPoint[2] = 0;
	SetPosition(0,0,0);
	SetFocalPoint(0,1,0);
	SetViewUp(0,0,1);
}

VtkTrackedCamera::~VtkTrackedCamera()
{
	QSettings settings("UFZ", "OpenGeoSys-5");
	settings.beginGroup("Tracking");
	settings.setValue("artOffsetX", _trackedPositionOffset[0]);
	settings.setValue("artOffsetY", _trackedPositionOffset[1]);
	settings.setValue("artOffsetZ", _trackedPositionOffset[2]);
	settings.setValue("artRealToVirtualScale", _realToVirtualScale);
	settings.setValue("artAspectRatio", _screenAspectRatio);
	settings.setValue("artScreenHeight", _screenHeight);
	settings.endGroup();
}

void VtkTrackedCamera::setTrackingData(double x, double y, double z)
{
	_trackedPosition[0] = x + _trackedPositionOffset[0];
	_trackedPosition[1] = y + _trackedPositionOffset[1];
	_trackedPosition[2] = z + _trackedPositionOffset[2];
	updateView();
}

void VtkTrackedCamera::setFocalPoint(double x, double y, double z)
{
	_focalPoint[0] = x;
	_focalPoint[1] = y;
	_focalPoint[2] = z;
	updateView();
}

void VtkTrackedCamera::setTrackingOffset(double x, double y, double z)
{
	setTrackingOffsetX(x);
	setTrackingOffsetY(y);
	setTrackingOffsetZ(z);
}

void VtkTrackedCamera::translate(double x, double y, double z)
{
	setFocalPoint(_focalPoint[0] + x, _focalPoint[1] + y, _focalPoint[2] + z);
	updateView();
}

void VtkTrackedCamera::rotate(double yaw, double elevation, double roll)
{
	Q_UNUSED(yaw);
	Q_UNUSED(elevation);
	Q_UNUSED(roll);
}

void VtkTrackedCamera::updateView()
{
	double x = _trackedPosition[0] * _realToVirtualScale;
	double y = _trackedPosition[1] * _realToVirtualScale;
	double z = _trackedPosition[2] * _realToVirtualScale;
	double angle = vtkMath::DegreesFromRadians(
	        2 * atan(0.5 * _screenHeight * _realToVirtualScale / y));

	//double newfocal[3] = {_x + x, _y, _z + z};
	//double newpos[3] = {_x + x, _y + y, _z + z};
	//double viewup[3] = {0, 0, 1};
	SetViewAngle(angle);
	//SetPosition(newpos);
	//SetFocalPoint(newfocal);
	//SetViewUp(viewup);
	//SetClippingRange(_zNear, _zFar);
	SetViewShear(x / y, z / y, 1); // see http://www.vtk.org/pipermail/vtkusers/2005-March/078735.html

	emit viewUpdated();
}

void VtkTrackedCamera::updatedFromOutside()
{
	double pos[3], dir[3];
	GetPosition(pos);
	GetDirectionOfProjection(dir);

	vtkMath::Normalize(dir);            // Get the view direction
	dir[0] = dir[0] * _trackedPosition[1]; // Multiplying the view direction
	dir[1] = dir[1] * _trackedPosition[1]; // with the tracked distance to the
	dir[2] = dir[2] * _trackedPosition[1]; // display results in the vector
	_focalPoint[0] = dir[0] + pos[0];   // from the position to the new
	_focalPoint[1] = dir[1] + pos[1];   // focal point.
	_focalPoint[2] = dir[2] + pos[2];

	updateView();

	//std::cout << "Camera slot" << std::endl;
}
