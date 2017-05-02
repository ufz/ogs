/**
 * \file
 * \author Lars Bilke
 * \date   2010-06-21
 * \brief  Definition of the VtkCustomInteractorStyle class.
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
#include <QObject>

#include <vtkInteractorStyleTrackballCamera.h>

class vtkDataObject;
class vtkDataSetMapper;
class vtkActor;

class vtkUnstructuredGridAlgorithm;

/**
 * VtkCustomInteractorStyle implements highlighting of an active actor and
 * highlighting of picked cells inside a vtk object.
 * Picking occurs when a vtk object was selected, the alternate mouse mode is
 * active (hold spacebar) and the user left clicks. On right click the cameras
 * focal point (center of rotation) is set to the picking position.
 */
class VtkCustomInteractorStyle : public QObject, public vtkInteractorStyleTrackballCamera
{
    Q_OBJECT

public:
    static VtkCustomInteractorStyle* New();
    vtkTypeMacro (VtkCustomInteractorStyle, vtkInteractorStyleTrackballCamera);

    /// @brief Handles key press events.
    void OnChar() override;

    /// @brief Handles key down events.
    void OnKeyDown() override;

    /// @brief Handles key up events.
    void OnKeyUp() override;

    /// @brief Handles left mouse button events (picking).
    void OnLeftButtonDown() override;

    /// @brief Handles middle mouse button events (rotation point picking).
    void OnRightButtonDown() override;

public slots:
    void highlightActor(vtkProp3D* prop);

    /// @brief Removes the highlight actor from the visible scene.
    void removeHighlightActor();

    void setHighlightActor(bool on);

    /// @brief Sets the highlightable vtk object.
    void pickableDataObject(vtkDataObject* object);

protected:
    VtkCustomInteractorStyle();
    ~VtkCustomInteractorStyle() override;

    /// @brief The vtk object to pick.
    vtkDataObject* _data;

    /// @brief The mapper for highlighting the selected cell.
    vtkDataSetMapper* _selectedMapper;

    /// @brief The actor for highlighting the selected cell.
    vtkActor* _selectedActor;

private:
    bool _highlightActor;
    bool _alternateMouseActions;

signals:
    /// @brief Emitted when something was picked.
    void requestViewUpdate();

    /// @brief Emitted when the cursor shape was changed due to alternate
    /// mouse action mode.
    void cursorChanged(Qt::CursorShape);

    /// @brief Emitted when a mesh element has been picked
    void elementPicked(vtkUnstructuredGridAlgorithm const*const, unsigned);

    /// @brief Emitted when the current object type cannot be handled by the element model
    void clearElementView();

};
