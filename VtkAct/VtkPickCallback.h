/**
 * \file VtkPickCallback.h
 * 21/6/2010 LB Initial implementation
 * 
 */


#ifndef VTKPICKCALLBACK_H
#define VTKPICKCALLBACK_H

// ** INCLUDES **
#include <QObject>
#include <vtkCommand.h>

class vtkProp3D;

/**
 * VtkPickCallback is a vtkCommand that implements functionality when
 * picking a vtk object through a vtkCellPicker.
 */
class VtkPickCallback : public QObject, public vtkCommand
{
	Q_OBJECT

public:
	static VtkPickCallback* New();

	void Execute(vtkObject *caller, unsigned long eventId, void *callData);	

protected:
	VtkPickCallback();

signals:
	/// Is emitted when an vtkActor was picked.
	void actorPicked (vtkProp3D* actor);

};

#endif // VTKPICKCALLBACK_H
