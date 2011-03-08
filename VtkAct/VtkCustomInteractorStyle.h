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
 */
class VtkCustomInteractorStyle : public QObject, public vtkInteractorStyleTrackballCamera
{
	Q_OBJECT

public:
	static VtkCustomInteractorStyle* New();
	vtkTypeMacro (VtkCustomInteractorStyle, vtkInteractorStyleTrackballCamera);

	/// @biref Handles key press events.
	virtual void OnChar();
	
	/// @brief Handles left mouse button events.
	virtual void OnLeftButtonDown();

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

};

#endif // VTKINTERACTORSTYLE_H
