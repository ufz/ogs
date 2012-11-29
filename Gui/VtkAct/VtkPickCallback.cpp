/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file VtkPickCallback.cpp
 *
 * Created on 2010-06-21 by Lars Bilke
 */

// ** INCLUDES **
#include "VtkPickCallback.h"

#include <vtkActor.h>
#include <vtkCellPicker.h>

VtkPickCallback* VtkPickCallback::New()
{
	return new VtkPickCallback();
}

void VtkPickCallback::Execute( vtkObject* caller, unsigned long vtkNotUsed(
                                       eventId), void* vtkNotUsed(callData) )
{
	vtkCellPicker* picker = static_cast<vtkCellPicker*>(caller);
	if (picker->GetCellId() < 0)
	{
		// Nothing is picked
	}
	else
	{
		vtkActor* actor = picker->GetActor();
		if (actor)
			emit actorPicked (actor);

		double* pos = picker->GetPickPosition();
		std::cout << "Picked cell id is: " << picker->GetCellId() << std::endl;
		std::cout << "Picked position is: " << pos[0] << " " << pos[1] << " " << pos[2] <<
		std::endl;
	}
}

VtkPickCallback::VtkPickCallback()
	: QObject()
{
}
