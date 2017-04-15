/**
 * \file
 * \author Lars Bilke
 * \date   2009-11-03
 * \brief  Definition of the VisualizationWidget class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

// ** INCLUDES **
#include "ui_VisualizationWidgetBase.h"

class vtkRenderer;
class vtkOrientationMarkerWidget;
class VtkCustomInteractorStyle;
class VtkPickCallback;

/**
 * \brief Widget containing the 3d VTK scene view.
 */
class VisualizationWidget : public QWidget, public Ui_VisualizationWidgetBase
{
    Q_OBJECT

public:

    /// @brief Constructor.
    VisualizationWidget(QWidget* parent = nullptr);

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
    /// x,y,z are in {-1, 0, 1} and specify from which direction the scene is displayed.
    void showAll(int x, int y, int z);

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

    /// @brief Resets the camera to view the entire scene.
    void on_showAllPushButton_pressed() { this->showAll(0,0,1); };

    /// @brief Reset camera to view entire scene from +x perspective.
    void on_rotateXPosPushButton_pressed() { this->showAll(1,0,0); };

    /// @brief Reset camera to view entire scene from -x perspective.
    void on_rotateXNegPushButton_pressed() { this->showAll(-1,0,0); };

    /// @brief Reset camera to view entire scene from +y perspective.
    void on_rotateYPosPushButton_pressed() { this->showAll(0,1,0); };

    /// @brief Reset camera to view entire scene from -y perspective.
    void on_rotateYNegPushButton_pressed() { this->showAll(0,-1,0); };

    /// @brief Reset camera to view entire scene from +z perspective.
    void on_rotateZPosPushButton_pressed() { this->showAll(0,0,1); };

    /// @brief Reset camera to view entire scene from -z perspective.
    void on_rotateZNegPushButton_pressed() { this->showAll(0,0,-1); };

    /// @brief Toggles rectangular zooming mode.
    void on_zoomToolButton_toggled(bool checked);

    /// @brief Toggles the display of bounding boxes around.
    void on_highlightToolButton_toggled(bool checked);

    /// @brief Toggles the orthogonal projection.
    void on_orthogonalProjectionToolButton_toggled(bool checked);

    /// @brief Saves a screenshot.
    void on_screenshotPushButton_pressed();

private:
    vtkRenderer* _vtkRender;
    vtkOrientationMarkerWidget* _markerWidget;
    VtkCustomInteractorStyle* _interactorStyle;
    VtkPickCallback* _vtkPickCallback;
    bool _isShowAllOnLoad;
};
