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

	/// Constructor.
	VisualizationWidget(QWidget* parent = 0);

	/// Destructor.
	~VisualizationWidget();

	/// Returns the VtkCustomInteractorStyle.
	VtkCustomInteractorStyle* interactorStyle() const;

	/// Returns the VtkPickCallback.
	VtkPickCallback* vtkPickCallback() const;

public slots:
	/// Updates the the 3d view.
	void updateView();

	/// Shows the entire scene on the views.
	void showAll();

	/// Returns the vtk renderer
	vtkRenderer* renderer() const { return _vtkRender; }

protected slots:
	/// Toggles stereo rendering on / off
	void on_stereoToolButton_toggled(bool checked);

	/// TODO Toggles full screen mode (actually not working)
	void on_fullscreenToolButton_clicked(bool checked);

	/// Adjusts the eye angle (separation) for stereo viewing
	void on_eyeAngleSlider_valueChanged(int value);

	/// Toggles rectangular zooming mode.
	void on_zoomToolButton_toggled(bool checked);
	
	/// Resets the camera to view the entire scene
	void on_showAllPushButton_pressed();
	
	/// Toggles the display of bounding boxes around
	void on_highlightToolButton_toggled(bool checked);

	/// Toggles the orthogonal projection
	void on_orthogonalProjectionToolButton_toggled(bool checked);

private:
	vtkRenderer* _vtkRender;
	VtkCustomInteractorStyle* _interactorStyle;
	VtkPickCallback* _vtkPickCallback;
	#ifdef OGS_USE_VRPN
	vtkEventQtSlotConnect* _qtConnect;
	#endif // OGS_USE_VRPN
};

#endif // VISUALIZATIONWIDGET_H
