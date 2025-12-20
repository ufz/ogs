// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

// ** INCLUDES **
#include "VtkPickCallback.h"

#include <vtkActor.h>
#include <vtkCellPicker.h>

#include "BaseLib/Logging.h"

VtkPickCallback* VtkPickCallback::New()
{
    return new VtkPickCallback();
}

void VtkPickCallback::Execute(vtkObject* caller,
                              unsigned long vtkNotUsed(eventId) /*eventId*/,
                              void* vtkNotUsed(callData) /*callData*/)
{
    auto* picker = static_cast<vtkCellPicker*>(caller);
    if (picker->GetCellId() < 0)
    {
        // Nothing is picked
    }
    else
    {
        vtkActor* actor = picker->GetActor();
        if (actor)
        {
            emit actorPicked(actor);
        }

        double* pos = picker->GetPickPosition();
        INFO("Picked cell id is: {:d}", picker->GetCellId());
        INFO("Picked position is: {:f} {:f} {:f}", pos[0], pos[1], pos[2]);
    }
}

VtkPickCallback::VtkPickCallback() {}
