/**
 * \file VtkApplyColorTableFilter.cpp
 * 21/10/2010 KR Initial implementation
 *
 */

// ** VTK INCLUDES **
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include "vtkObjectFactory.h"
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkLookupTable.h>
#include <vtkPointData.h>
#include <vtkCellData.h>

#include "VtkApplyColorTableFilter.h"

vtkStandardNewMacro(VtkApplyColorTableFilter);
vtkCxxSetObjectMacro(VtkApplyColorTableFilter, ColorLookupTable, vtkLookupTable);
vtkCxxRevisionMacro(VtkApplyColorTableFilter, "$Revision: 6575 $");


VtkApplyColorTableFilter::VtkApplyColorTableFilter() : ColorLookupTable(NULL)
{
	this->SetColorsOnCells(false);
}

VtkApplyColorTableFilter::~VtkApplyColorTableFilter()
{
}

int VtkApplyColorTableFilter::RequestData( vtkInformation* request, 
							             vtkInformationVector** inputVector, 
								         vtkInformationVector* outputVector )
{
	if (this->ColorLookupTable==NULL) return 0;

	(void)request;

	vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
    vtkPolyData *input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

	vtkSmartPointer<vtkUnsignedCharArray> colorTable = this->ColorLookupTable->GetTable();
	vtkSmartPointer<vtkUnsignedCharArray> colorArray = vtkSmartPointer<vtkUnsignedCharArray>::New();
		colorArray->SetNumberOfComponents(4);
		colorArray->SetName("Colors");
	vtkSmartPointer<vtkUnsignedCharArray> scalars = 
		(!ColorsOnCells) ? vtkUnsignedCharArray::SafeDownCast(input->GetPointData()->GetScalars())
						  : vtkUnsignedCharArray::SafeDownCast(input->GetCellData()->GetScalars());
	int limit = (!ColorsOnCells) ? input->GetNumberOfPoints() : input->GetNumberOfCells();
	
	for (int i=0; i<limit; i++)
	{
		double* value = scalars->GetTuple(i);
		unsigned char* rgba = new unsigned char[4];
		colorTable->GetTupleValue(static_cast<int>(value[0]), rgba);
		colorArray->InsertNextTupleValue(rgba); 
	}

	vtkInformation* outInfo = outputVector->GetInformationObject(0);
    vtkPolyData* output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
		output->CopyStructure(input);
		output->GetPointData()->PassData(input->GetPointData());
		output->GetCellData() ->PassData(input->GetCellData());
	if (!ColorsOnCells)
	{
		output->GetPointData()->AddArray(colorArray);
		output->GetPointData()->SetActiveScalars("Colors");
	}
	else
	{
		output->GetCellData()->AddArray(colorArray);
		output->GetCellData()->SetActiveScalars("Colors");
	}

	this->ColorLookupTable->Build();

	return 1;
}

unsigned long VtkApplyColorTableFilter::GetMTime()
{
	unsigned long t1, t2;

	t1 = this->Superclass::GetMTime();
	if (this->ColorLookupTable)
	{
		t2 = this->ColorLookupTable->GetMTime();
		if (t2 > t1)
			t1 = t2;
	}
	return t1;
}

void VtkApplyColorTableFilter::PrintSelf( ostream& os, vtkIndent indent )
{
	this->Superclass::PrintSelf(os, indent);
	os << indent << "ColorsOnCells: " << this->ColorsOnCells << "\n";
	ColorLookupTable->PrintSelf(os, indent);
}
