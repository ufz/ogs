/**
 * \file VtkColorByHeightFilter.cpp
 * 21/04/2010 KR Initial implementation
 *
 */

// ** VTK INCLUDES **
#include "VtkColorByHeightFilter.h"
#include "VtkColorLookupTable.h"

#include <vtkCellData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkLookupTable.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkStreamingDemandDrivenPipeline.h>

vtkStandardNewMacro(VtkColorByHeightFilter);
vtkCxxRevisionMacro(VtkColorByHeightFilter, "$Revision$");

VtkColorByHeightFilter::VtkColorByHeightFilter()
{
	ColorLookupTable = VtkColorLookupTable::New();
	ColorLookupTable->GetTableRange(this->_tableRange);
	ColorLookupTable->setInterpolationType(VtkColorLookupTable::LINEAR);
	_tableRangeScaling = 1.0;
	_activeAttributeName = "P-Colors";
}

VtkColorByHeightFilter::~VtkColorByHeightFilter()
{
	ColorLookupTable->Delete();
}

void VtkColorByHeightFilter::PrintSelf( ostream& os, vtkIndent indent )
{
	this->Superclass::PrintSelf(os,indent);

	double range[2];
	ColorLookupTable->GetTableRange(range);
	os << indent << "== VtkColorByHeightFilter ==" << endl;
	os << indent << "Range: " << range[0] << "-" << range[1] << endl;
	os << indent << "Interpolation Type:" << ColorLookupTable->getInterpolationType() << endl;
}

unsigned long VtkColorByHeightFilter::GetMTime()
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
int VtkColorByHeightFilter::RequestData( vtkInformation*,
                                         vtkInformationVector** inputVector,
                                         vtkInformationVector* outputVector )
{
	vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
	vtkPolyData* input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

	ColorLookupTable->SetTableRange(getMinHeight(input), getMaxHeight(input));
	ColorLookupTable->Build();

	vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->SetName("Colors");

	// Generate the colors for each point based on the color map
	size_t nPoints = input->GetNumberOfPoints();
	for (size_t i = 0; i < nPoints; i++)
	{
		double p[3];
		input->GetPoint(i,p);

		unsigned char lutColor[4];
		ColorLookupTable->getColor((int)p[2], lutColor);
		colors->InsertNextTupleValue(lutColor);
	}

	vtkInformation* outInfo = outputVector->GetInformationObject(0);
	vtkPolyData* output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));
	output->CopyStructure(input);
	output->GetPointData()->PassData(input->GetPointData());
	output->GetCellData()->PassData(input->GetCellData());
	output->GetPointData()->AddArray(colors);
	output->GetPointData()->SetActiveScalars("Colors");

	return 1;
}

double VtkColorByHeightFilter::getMinHeight(vtkPolyData* data)
{
	double range[2];
	ColorLookupTable->GetTableRange(range);
	double min = range[0];
	size_t nPoints = data->GetNumberOfPoints();
	if (min == ColorLookupTable->DEFAULTMINVALUE && nPoints > 0)
	{
		double p[3];
		data->GetPoint(0,p);
		min = p[2];

		for (size_t i = 1; i < nPoints; i++)
		{
			data->GetPoint(i,p);
			if (p[2] < min)
				min = p[2];
		}
	}
	return min;
}

double VtkColorByHeightFilter::getMaxHeight(vtkPolyData* data)
{
	double range[2];
	ColorLookupTable->GetTableRange(range);
	double max = range[1];
	size_t nPoints = data->GetNumberOfPoints();
	if (max == ColorLookupTable->DEFAULTMAXVALUE && nPoints > 0)
	{
		double p[3];
		data->GetPoint(0,p);
		max = p[2];

		for (size_t i = 1; i < nPoints; i++)
		{
			data->GetPoint(i,p);
			if (p[2] > max)
				max = p[2];
		}
	}
	return max;
}

void VtkColorByHeightFilter::SetTableRange(double min, double max)
{
	if (min < max)
	{
		this->_tableRange[0] = min;
		this->_tableRange[1] = max;
		this->ColorLookupTable->SetTableRange(min, max);
	}
	else
		vtkstd::cout <<
		"VtkColorByHeightFilter::SetLimits(min, max) - Limits not changed because min value > max value."
		             << vtkstd::endl;
}

void VtkColorByHeightFilter::SetTableRangeScaling( double scale )
{
	this->_tableRangeScaling = scale;
	this->ColorLookupTable->SetTableRange(
	        this->_tableRange[0] * scale, this->_tableRange[1] * scale);
}
