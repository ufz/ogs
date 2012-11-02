/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file VisualizationWidget.h
 *
 * Created on 2009-11-03 by Lars Bilke
 */

#ifndef VISUALIZATIONWIDGET_H
#define VISUALIZATIONWIDGET_H

// ** INCLUDES **
#include "ui_VisualizationWidgetBase.h"

class vtkRenderer;
class VtkCustomInteractorStyle;
class VtkPickCallback;
#ifdef OGS_USE_VRPN
class vtkEventQtSlotConnect;
#endif // OGS_USE_VRPN

/**
 * \brief Widget containing the 3d VTK scene view.
 */
class VisualizationWidget : public QWidget, public Ui_VisualizationWidgetBase
{
	Q_OBJECT

public:

	/// @brief Constructor.
	VisualizationWidget(QWidget* parent = 0);

	/// @brief Destructor.
	~VisualizationWidget();

	/// @brief Returns the VtkCustomInteractorStyle.
	VtkCustomInteractorStyle* interactorStyle() const;

	/// @brief Returns the VtkPickCallback.
	VtkPickCallback* vtkPickCallback() const;

	/// @brief See updateViewOnLoad().
	void setShowAllOnLoad(bool show) { _isShowAllOnLoad = show; }

public slots:
	/// @brief Updates the the 3d view.
	void updateView();

	/// @brief Shows the entire scene on the views.
	void showAll();

	/// @brief Updates the view only or additionally shows the entire scene.
	void updateViewOnLoad();

	/// @brief Saves a magnified image of the current render window to a file.
	void screenshot(QString filename, int magnification);

	/// @brief Returns the vtk renderer
	vtkRenderer* renderer() const { return _vtkRender; }

	/// @brief Sets the widgets cursor shape.
	/// @see http://doc.qt.nokia.com/4.7/qt.html#CursorShape-enum
	void setCursorShape(Qt::CursorShape shape);

protected slots:

	/// @brief Toggles rectangular zooming mode.
	void on_zoomToolButton_toggled(bool checked);

	/// @brief Resets the camera to view the entire scene.
	void on_showAllPushButton_pressed();

	/// @brief Toggles the display of bounding boxes around.
	void on_highlightToolButton_toggled(bool checked);

	/// @brief Toggles the orthogonal projection.
	void on_orthogonalProjectionToolButton_toggled(bool checked);

	/// @brief Saves a screenshot.
	void on_screenshotPushButton_pressed();

private:
	vtkRenderer* _vtkRender;
	VtkCustomInteractorStyle* _interactorStyle;
	VtkPickCallback* _vtkPickCallback;
	bool _isShowAllOnLoad;
#ifdef OGS_USE_VRPN
	vtkEventQtSlotConnect* _qtConnect;
#endif     // OGS_USE_VRPN
};

#endif // VISUALIZATIONWIDGET_H
