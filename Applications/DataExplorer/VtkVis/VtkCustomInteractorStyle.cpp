/**
 * \file
 * \author Lars Bilke
 * \date   2010-06-21
 * \brief  Implementation of the VtkCustomInteractorStyle class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkCustomInteractorStyle.h"

#include <logog/include/logog.hpp>

#include <vtkActor.h>
#include <vtkAlgorithmOutput.h>
#include <vtkCamera.h>
#include <vtkCellData.h>
#include <vtkCellPicker.h>
#include <vtkDataSetMapper.h>
#include <vtkExtractSelection.h>
#include <vtkIdTypeArray.h>
#include <vtkObjectFactory.h>
#include <vtkProp.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRendererCollection.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridAlgorithm.h>

#include <string>

#include "VtkCompositeElementSelectionFilter.h"

vtkStandardNewMacro(VtkCustomInteractorStyle);

VtkCustomInteractorStyle::VtkCustomInteractorStyle()
: _data(nullptr), _highlightActor(false), _alternateMouseActions(false)
{
    _selectedMapper = vtkDataSetMapper::New();
    _selectedActor = vtkActor::New();
    _selectedActor->SetMapper(_selectedMapper);
    _selectedActor->GetProperty()->EdgeVisibilityOn();
    _selectedActor->GetProperty()->SetEdgeColor(1,0,0);
    _selectedActor->GetProperty()->SetLineWidth(3);
}

VtkCustomInteractorStyle::~VtkCustomInteractorStyle()
{
    _selectedActor->Delete();
    _selectedMapper->Delete();
}

void VtkCustomInteractorStyle::OnChar()
{
    switch (Interactor->GetKeyCode())
    {
    case '3':
        INFO("The 3 key was pressed.");
        break;
    case 'a':
        break;
    default:
        vtkInteractorStyleTrackballCamera::OnChar();
    }
}

void VtkCustomInteractorStyle::OnKeyDown()
{
    switch (Interactor->GetKeyCode())
    {
    case 32: // Space
        _alternateMouseActions = true;
        emit cursorChanged(Qt::CrossCursor);
        break;
    default:
        vtkInteractorStyleTrackballCamera::OnKeyDown();
    }
}

void VtkCustomInteractorStyle::OnKeyUp()
{
    switch (Interactor->GetKeyCode())
    {
    case 32: // Space
        _alternateMouseActions = false;
        emit cursorChanged(Qt::ArrowCursor);
        break;
    default:
        vtkInteractorStyleTrackballCamera::OnKeyUp();
    }
}

void VtkCustomInteractorStyle::highlightActor( vtkProp3D* actor )
{
    if (_highlightActor)
        HighlightProp((vtkProp*)actor);
}

void VtkCustomInteractorStyle::removeHighlightActor()
{
    this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->
        RemoveActor(_selectedActor);
}

void VtkCustomInteractorStyle::setHighlightActor(bool on)
{
    _highlightActor = on;
    if (!on)
        HighlightProp((vtkProp*)nullptr);
}

void VtkCustomInteractorStyle::pickableDataObject(vtkDataObject* object)
{
    _data = object;
    if (!object)
    {
        this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->
        RemoveActor(_selectedActor);
        _selectedMapper->SetInputConnection(nullptr);
    }
}

// From http://www.vtk.org/Wiki/VTK/Examples/Cxx/Picking/CellPicking
void VtkCustomInteractorStyle::OnLeftButtonDown()
{
    if (!_data)
        return vtkInteractorStyleTrackballCamera::OnLeftButtonDown();

    if (_alternateMouseActions)
    {
        // Get the location of the click (in window coordinates)
        int* pos = this->GetInteractor()->GetEventPosition();

        vtkSmartPointer<vtkCellPicker> picker =
                vtkSmartPointer<vtkCellPicker>::New();
        picker->SetTolerance(0.0005);

        // Pick from this location.
        picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());

        double* worldPosition = picker->GetPickPosition();
        INFO("Cell id is: %d", picker->GetCellId());

        if(picker->GetCellId() != -1)
        {
            INFO("Pick position is: %f %f %f", worldPosition[0], worldPosition[1], worldPosition[2]);

            vtkSmartPointer<vtkIdTypeArray> ids =
                    vtkSmartPointer<vtkIdTypeArray>::New();
            ids->SetNumberOfValues(1);
            ids->SetValue(0, picker->GetCellId());

            vtkSmartPointer<vtkSelectionNode> selectionNode =
                    vtkSmartPointer<vtkSelectionNode>::New();
            selectionNode->SetFieldType(vtkSelectionNode::CELL);
            selectionNode->SetContentType(vtkSelectionNode::INDICES);
            selectionNode->SetSelectionList(ids);

            vtkSmartPointer<vtkSelection> selection =
                    vtkSmartPointer<vtkSelection>::New();
            selection->AddNode(selectionNode);

            vtkSmartPointer<vtkExtractSelection> extractSelection =
                    vtkSmartPointer<vtkExtractSelection>::New();
            extractSelection->SetInputData(0, _data);
            extractSelection->SetInputData(1, selection);
            extractSelection->Update();

            // In selection
            vtkSmartPointer<vtkUnstructuredGrid> selected =
                    vtkSmartPointer<vtkUnstructuredGrid>::New();
            selected->ShallowCopy(extractSelection->GetOutput());

            INFO("There are %d points in the selection.", selected->GetNumberOfPoints());
            INFO("There are %d cells in the selection.", selected->GetNumberOfCells());

            // check if the underlying object is a mesh and if so, send a signal to the element model for display of information about the picked element.
            vtkAlgorithm* data_set = picker->GetActor()->GetMapper()->GetInputConnection(0, 0)->GetProducer()->GetInputConnection(0,0)->GetProducer();
            auto* source =
                dynamic_cast<vtkUnstructuredGridAlgorithm*>(data_set);
            if (source)
                emit elementPicked(source, static_cast<unsigned>(picker->GetCellId()));
            else
                emit clearElementView();
            _selectedMapper->SetInputData(selected);

            this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->
            AddActor(_selectedActor);
            //_highlightActor = true;
        }
        else
            this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->
            RemoveActor(_selectedActor);
        emit requestViewUpdate();
    }
    else
        // Forward events
        vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
}

void VtkCustomInteractorStyle::OnRightButtonDown()
{
    if (!_data)
        return vtkInteractorStyleTrackballCamera::OnRightButtonDown();

    if (_alternateMouseActions)
    {
        // Get the location of the click (in window coordinates)
        int* pos = this->GetInteractor()->GetEventPosition();

        vtkSmartPointer<vtkCellPicker> picker =
                vtkSmartPointer<vtkCellPicker>::New();
        picker->SetTolerance(0.0005);

        // Pick from this location.
        picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());

        double* worldPosition = picker->GetPickPosition();
        INFO("Cell id is: %d", picker->GetCellId());

        if(picker->GetCellId() != -1)
        {
            vtkRenderer* renderer =
                    this->Interactor->GetRenderWindow()->GetRenderers()->
                    GetFirstRenderer();
            vtkCamera* cam = renderer->GetActiveCamera();
            cam->SetFocalPoint(worldPosition);
            emit requestViewUpdate();
        }
    }
    else
        // Forward events
        vtkInteractorStyleTrackballCamera::OnRightButtonDown();
}
