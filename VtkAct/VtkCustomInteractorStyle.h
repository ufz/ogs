/**
 * \file VtkCustomInteractorStyle.h
 * 21/6/2010 LB Initial implementation
 * 
 */


#ifndef VTKCUSTOMINTERACTORSTYLE_H
#define VTKCUSTOMINTERACTORSTYLE_H

// ** INCLUDES **
#include <QObject>

#include <vtkInteractorStyleTrackballCamera.h>

class vtkDataObject;
class vtkDataSetMapper;
class vtkActor;

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

	/// @biref Handles key press events.
	virtual void OnChar();

	/// @biref Handles key down events.
	virtual void OnKeyDown();

	/// @biref Handles key up events.
	virtual void OnKeyUp();
	
	/// @brief Handles left mouse button events (picking).
	virtual void OnLeftButtonDown();

	/// @brief Handles middle mouse button events (rotation point picking).
	virtual void OnRightButtonDown();

public slots:
	void highlightActor(vtkProp3D* prop);
	void setHighlightActor(bool on);
	
	/// @brief Sets the highlightable vtk object.
	void pickableDataObject(vtkDataObject* object);

protected:
	VtkCustomInteractorStyle();
	virtual ~VtkCustomInteractorStyle();
	
	/// @brief The vtk object to pick.
	vtkDataObject* Data;
	
	/// @brief The mapper for highlighting the selected cell.
	vtkDataSetMapper* selectedMapper;
	
	/// @brief The actor for highlighting the selected cell.
	vtkActor* selectedActor;
	

private:
	bool _highlightActor;
	bool _alternateMouseActions;

signals:
	/// @brief Emitted when something was picked.
	void requestViewUpdate();

	/// @brief Emitted when the cursor shape was changed due to alternate
	/// mouse action mode.
	void cursorChanged(Qt::CursorShape);

};

#endif // VTKINTERACTORSTYLE_H
