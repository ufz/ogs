/**
 * \file VtkConditionSource.cpp
 * 2011/03/02 KR Initial implementation
 */

// ** INCLUDES **
#include "VtkConditionSource.h"

#include <vtkCellArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkSmartPointer.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolygon.h>

vtkStandardNewMacro(VtkConditionSource);
vtkCxxRevisionMacro(VtkConditionSource, "$Revision$");

VtkConditionSource::VtkConditionSource()
: _points(NULL), _polylines(NULL), _surfaces(NULL)
{
	this->SetNumberOfInputPorts(0);

	//const GEOLIB::Color* c = GEOLIB::getRandomColor();
	//GetProperties()->SetColor((*c)[0]/255.0,(*c)[1]/255.0,(*c)[2]/255.0);
}

void VtkConditionSource::setData(const std::vector<GEOLIB::Point*>* points, const std::vector<GEOLIB::Polyline*>* lines, const std::vector<GEOLIB::Surface*>* surfaces)
{
	_points    = points;
	_polylines = lines;
	_surfaces  = surfaces;
}

void VtkConditionSource::PrintSelf( ostream& os, vtkIndent indent )
{
	this->Superclass::PrintSelf(os,indent);

	if (_points->size() == 0)
		return;

	os << indent << "== VtkConditionSource ==" << "\n";

}

int VtkConditionSource::RequestData( vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector )
{
	(void)request;
	(void)inputVector;

	if 	(_points ) 
	{
		if (_points->empty())
		{
			std::cout << "ERROR in VtkConditionSource::RequestData : Size of point vector is 0" << std::endl;
			return 0;
		}
	}
	else return 0;

	vtkSmartPointer<vtkInformation> outInfo = outputVector->GetInformationObject(0);
	vtkSmartPointer<vtkPolyData> output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

	vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> newCells = vtkSmartPointer<vtkCellArray>::New();

	if (outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()) > 0)
		return 1;

	vtkSmartPointer<vtkIntArray> plyIDs = vtkSmartPointer<vtkIntArray>::New();
		plyIDs->SetNumberOfComponents(1);
		plyIDs->SetName("PolylineIDs");
	vtkSmartPointer<vtkIntArray> sfcIDs = vtkSmartPointer<vtkIntArray>::New();
		sfcIDs->SetNumberOfComponents(1);
		sfcIDs->SetName("SurfaceIDs");

	
	for (std::vector<GEOLIB::Point*>::const_iterator it = _points->begin();
		it != _points->end(); ++it)
	{
		double coords[3] = {(*(*it))[0], (*(*it))[1], (*(*it))[2]};
		vtkIdType pid = newPoints->InsertNextPoint(coords);
		newCells->InsertNextCell(1, &pid);
	}

	if (!_polylines->empty())
	{
		//int lastMaxIndex = 0;
		for (size_t j=0; j<_polylines->size(); j++)
		{
			const int nPoints = (*_polylines)[j]->getNumberOfPoints();
/*
			for (int i = 0; i < numPoints; i++)
			{
				const double* coords = (*(*_polylines)[j])[i]->getData();
				newPoints->InsertNextPoint(coords);
			}
*/
			// Generate lines
			newCells->InsertNextCell(nPoints);
			//plyIDs->InsertNextValue(j);
			for (int i = 0; i < nPoints; i++)
				newCells->InsertCellPoint((*(*_polylines)[j]).getPointID(i));

			//lastMaxIndex += numPoints;
		}
	}

	if (!_surfaces->empty())
	{
/*
		for (size_t i=0; i<nPoints; i++)
		{
			double* coords = const_cast<double*>((*surfacePoints)[i]->getData());
			newPoints->InsertNextPoint(coords);
		}
*/
		//vtkIdType count(0);
		for (std::vector<GEOLIB::Surface*>::const_iterator it = _surfaces->begin();
			it != _surfaces->end(); ++it)
		{
			const size_t nTriangles = (*it)->getNTriangles();

			for (size_t i = 0; i < nTriangles; i++)
			{
				vtkPolygon* aPolygon = vtkPolygon::New();
				aPolygon->GetPointIds()->SetNumberOfIds(3);

				const GEOLIB::Triangle* triangle = (**it)[i];
				for (size_t j=0; j<3; j++)
				{
					aPolygon->GetPointIds()->SetId(j, ((*triangle)[j]));
				}
				newCells->InsertNextCell(aPolygon);
				//sfcIDs->InsertNextValue(count);

				aPolygon->Delete();
			}
			//count++;
		}
	}
/*
	output->SetPoints(newPoints);
	output->SetPolys(newPolygons);
	output->GetCellData()->AddArray(sfcIDs);
	output->GetCellData()->SetActiveAttribute("SurfaceIDs", vtkDataSetAttributes::SCALARS);



	output->SetPoints(newPoints);
	output->SetLines(newLines);
	output->GetCellData()->AddArray(plyIDs);
	output->GetCellData()->SetActiveAttribute("PolylineIDs", vtkDataSetAttributes::SCALARS);
*/
	output->SetPoints(newPoints);
	output->SetVerts(newCells);

	return 1;
}

int VtkConditionSource::RequestInformation( vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector )
{
	(void)request;
	(void)inputVector;

	vtkInformation* outInfo = outputVector->GetInformationObject(0);
	outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);

	return 1;
}

void VtkConditionSource::SetUserProperty( QString name, QVariant value )
{
	Q_UNUSED(name);
	Q_UNUSED(value);
}
