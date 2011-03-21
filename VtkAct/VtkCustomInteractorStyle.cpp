/**
 * \file VtkInteractorStyle.cpp
 * 21/6/2010 LB Initial implementation
 * 
 * Implementation of VtkInteractorStyle
 */

// ** INCLUDES **
#include "VtkCustomInteractorStyle.h"

#include <vtkRenderWindowInteractor.h>
#include <vtkObjectFactory.h>
#include <vtkProp.h>
#include <vtkSmartPointer.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkCellPicker.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkExtractSelection.h>
#include <vtkIdTypeArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRendererCollection.h>

#include <string>

vtkStandardNewMacro(VtkCustomInteractorStyle);

VtkCustomInteractorStyle::VtkCustomInteractorStyle()
: _highlightActor(true)
{
	selectedMapper = vtkDataSetMapper::New();
	selectedActor = vtkActor::New();
	selectedActor->SetMapper(selectedMapper);
	selectedActor->GetProperty()->EdgeVisibilityOn();
	selectedActor->GetProperty()->SetEdgeColor(1,0,0);
	selectedActor->GetProperty()->SetLineWidth(3);
	Data = NULL;
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
		std::cout << "The 3 key was pressed." << std::endl;
		break;
	case 'a':
		break;
	default:
		vtkInteractorStyleTrackballCamera::OnChar();
	}
}

void VtkCustomInteractorStyle::highlightActor( vtkProp3D* actor )
{
	if (_highlightActor)
		HighlightProp((vtkProp*)actor);
}

void VtkCustomInteractorStyle::setHighlightActor(bool on)
{ 
	_highlightActor = on;
	if (!on)
		HighlightProp((vtkProp*)NULL);
}

void VtkCustomInteractorStyle::pickableDataObject(vtkDataObject* object)
{
	Data = object;
	if (!object)
	{
		this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(selectedActor);
		selectedMapper->SetInputConnection(NULL);
	}
}

// From http://www.vtk.org/Wiki/VTK/Examples/Cxx/Picking/CellPicking
void VtkCustomInteractorStyle::OnLeftButtonDown()
{
	if (!Data)
		return vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
	
	// Get the location of the click (in window coordinates)
	int* pos = this->GetInteractor()->GetEventPosition();
    
	vtkSmartPointer<vtkCellPicker> picker =
	  vtkSmartPointer<vtkCellPicker>::New();
	picker->SetTolerance(0.0005);
    
	// Pick from this location.
	picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());
    
	double* worldPosition = picker->GetPickPosition();
	std::cout << "Cell id is: " << picker->GetCellId() << std::endl;
    
	if(picker->GetCellId() != -1)
	  {
    
	  std::cout << "Pick position is: " << worldPosition[0] << " " << worldPosition[1]
	            << " " << worldPosition[2] << endl;
    
	  vtkSmartPointer<vtkIdTypeArray> ids =
	    vtkSmartPointer<vtkIdTypeArray>::New();
	  ids->SetNumberOfComponents(1);
	  ids->InsertNextValue(picker->GetCellId());
    
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
	  extractSelection->SetInput(0, this->Data);
	  extractSelection->SetInput(1, selection);
	  extractSelection->Update();
    
	  // In selection
	  vtkSmartPointer<vtkUnstructuredGrid> selected =
	    vtkSmartPointer<vtkUnstructuredGrid>::New();
	  selected->ShallowCopy(extractSelection->GetOutput());
    
	  std::cout << "There are " << selected->GetNumberOfPoints()
	            << " points in the selection." << std::endl;
	  std::cout << "There are " << selected->GetNumberOfCells()
	            << " cells in the selection." << std::endl;
    
    
	  selectedMapper->SetInputConnection(selected->GetProducerPort());
    
	  this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(selectedActor);
    
	  }
	// Forward events
	vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
}
