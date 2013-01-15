/**
 * \file
 * \author Lars Bilke
 * \date   2010-06-21
 * \brief  Implementation of the VtkPickCallback class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
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
