/**
 * \file
 * \author Lars Bilke
 * \date   2010-06-21
 * \brief  Implementation of the VtkCustomInteractorStyle class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "VtkCustomInteractorStyle.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

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
	: _highlightActor(true), _alternateMouseActions(false)
{
	selectedMapper = vtkDataSetMapper::New();
	selectedActor = vtkActor::New();
	selectedActor->SetMapper(selectedMapper);
	selectedActor->GetProperty()->EdgeVisibilityOn();
	selectedActor->GetProperty()->SetEdgeColor(1,0,0);
	selectedActor->GetProperty()->SetLineWidth(3);
	Data = nullptr;
}

VtkCustomInteractorStyle::~VtkCustomInteractorStyle()
{
	selectedActor->Delete();
	selectedMapper->Delete();
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
		RemoveActor(selectedActor);
}

void VtkCustomInteractorStyle::setHighlightActor(bool on)
{
	_highlightActor = on;
	if (!on)
		HighlightProp((vtkProp*)nullptr);
}

void VtkCustomInteractorStyle::pickableDataObject(vtkDataObject* object)
{
	Data = object;
	if (!object)
	{
		this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->
		RemoveActor(selectedActor);
		selectedMapper->SetInputConnection(nullptr);
	}
}

// From http://www.vtk.org/Wiki/VTK/Examples/Cxx/Picking/CellPicking
void VtkCustomInteractorStyle::OnLeftButtonDown()
{
	if (!Data)
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
			extractSelection->SetInputData(0, this->Data);
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
			vtkUnstructuredGridAlgorithm* source = dynamic_cast<vtkUnstructuredGridAlgorithm*>(data_set);
			if (source)
				emit elementPicked(source, static_cast<unsigned>(picker->GetCellId()));
			else
				emit clearElementView();
			selectedMapper->SetInputData(selected);

			this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->
			AddActor(selectedActor);
			_highlightActor = true;
		}
		else
			this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->
			RemoveActor(selectedActor);
		emit requestViewUpdate();
	}
	else
		// Forward events
		vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
}

void VtkCustomInteractorStyle::OnRightButtonDown()
{
	if (!Data)
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
