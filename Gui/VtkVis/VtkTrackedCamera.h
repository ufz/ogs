/**
 * \file
 * \author Lars Bilke
 * \date   2010-08-25
 * \brief  Definition of the VtkTrackedCamera class.
 *
 * \copyright
 * Copyright (c)  2013, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef VTKTRACKEDCAMERA_H
#define VTKTRACKEDCAMERA_H

#include <QObject>
#include <vtkOpenGLCamera.h>

/// @brief This camera implements view shearing for using with a head tracking
/// system.
class VtkTrackedCamera : public QObject, public vtkOpenGLCamera
{
	Q_OBJECT

public:
	/// @brief Constructor.
	VtkTrackedCamera(QObject* parent);

	/// @brief Destructor
	virtual ~VtkTrackedCamera();

	/// Sets the scaling from real meters to virtual meters.
	void setRealToVirtualScale(double scale) { _realToVirtualScale = scale; }
	double realToVirtualScale () const { return _realToVirtualScale; }

	double trackingOffsetX() const { return _trackedPositionOffset[0]; }
	double trackingOffsetY() const { return _trackedPositionOffset[1]; }
	double trackingOffsetZ() const { return _trackedPositionOffset[2]; }

	double screenAspectRatio() const { return _screenAspectRatio; }
	double screenHeight() const { return _screenHeight; }

public slots:
	//void setTrackinData(double position[3], double dir[3]);
	/// @brief Sets the tracked position. The coordinate origin must be the center
	/// of the projection wall. See also setTrackingOffset().
	void setTrackingData(double x, double y, double z);

	/// Sets the focal point in world coordinates. This point corresponds to the
	/// center of the projection wall.
	void setFocalPoint(double x, double y, double z);

	/// @brief Sets an offset of the tracked position. This must be used to setup
	/// the origin of the calibrated tracking space.
	void setTrackingOffset(double x, double y, double z);

	void setTrackingOffsetX(double val) { _trackedPositionOffset[0] = val; }
	void setTrackingOffsetY(double val) { _trackedPositionOffset[1] = val; }
	void setTrackingOffsetZ(double val) { _trackedPositionOffset[2] = val; }

	/// Move the camera by the given vector.
	void translate(double x, double y, double z);

	/// Rotate the camera by the given angles.
	void rotate(double yaw, double elevation, double roll);

	/// Updates the view.
	void updateView();

	/// Must be called to update the view after the camera was modified from the
	/// outside, e.g. from the vtkRenderWindowInteractor.
	void updatedFromOutside();

	/// @brief Sets the screen properties.
	/// @param aspectRatio width / height
	/// @param height The screen height in meter.
	void setScreenProperties(double aspectRatio, double height)
	{
		_screenAspectRatio = aspectRatio;
		_screenHeight = height;
	}

	/// @brief Sets the screen aspect ratio (width / height).
	void setScreenAspectRatio(double ratio) { _screenAspectRatio = ratio;
		                                  updateView(); }

	/// @brief Sets the screen height in meter.
	void setScreenHeight(double height) { _screenHeight = height;
		                              updateView(); }

private:
	double _focalPoint[3];
	double _trackedPosition[3];
	double _trackedPositionOffset[3];
	double _realToVirtualScale;
	double _screenAspectRatio;
	double _screenHeight;

signals:
	/// Is emitted from updateView().
	void viewUpdated();
};

#endif // VTKTRACKEDCAMERA_H
