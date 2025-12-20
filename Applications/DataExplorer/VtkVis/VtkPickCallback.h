// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

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

    void Execute(vtkObject* caller, unsigned long eventId,
                 void* callData) override;

protected:
    VtkPickCallback();

signals:
    /// Is emitted when an vtkActor was picked.
    void actorPicked (vtkProp3D* actor);
};
