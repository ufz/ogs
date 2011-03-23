	/**
 * \file VisualizationWidget.h
 * 3/11/2009 LB Initial implementation
 *
 */


#ifndef VISUALIZATIONWIDGET_H
#define VISUALIZATIONWIDGET_H

// ** INCLUDES **
#include "ui_VisualizationWidgetBase.h"
#include "Configure.h"

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

protected slots:
	/// @brief Toggles stereo rendering on / off.
	void on_stereoToolButton_toggled(bool checked);

	/// @brief Adjusts the eye angle (separation) for stereo viewing.
	//void on_eyeAngleSlider_valueChanged(int value);

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
	#endif // OGS_USE_VRPN
};

#endif // VISUALIZATIONWIDGET_H
