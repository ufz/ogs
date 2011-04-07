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

void VtkConditionSource::setData(const std::vector<GEOLIB::Point*>* points, const std::vector<GEOLIB::Polyline*>* lines, const std::vector<GEOLIB::Surface*>* surfaces,
								 std::vector<size_t> *pnt_idx, std::vector<size_t> *ply_idx, std::vector<size_t> *sfc_idx)
{
	_points    = points;
	_polylines = lines;
	_surfaces  = surfaces;
	_pnt_idx   = pnt_idx;
	_ply_idx   = ply_idx;
	_sfc_idx   = sfc_idx;
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

	if 	( _points ) 
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

	vtkSmartPointer<vtkPoints> newPoints   = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> newVerts = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkCellArray> newLines = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkCellArray> newPolys = vtkSmartPointer<vtkCellArray>::New();

	if (outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()) > 0)
		return 1;

	vtkSmartPointer<vtkIntArray> plyIDs = vtkSmartPointer<vtkIntArray>::New();
		plyIDs->SetNumberOfComponents(1);
		plyIDs->SetName("PolylineIDs");
	vtkSmartPointer<vtkIntArray> sfcIDs = vtkSmartPointer<vtkIntArray>::New();
		sfcIDs->SetNumberOfComponents(1);
		sfcIDs->SetName("SurfaceIDs");

	size_t nPoints = _points->size();
	for (size_t i=0; i<nPoints; i++)
	{
		double coords[3] = {(*(*_points)[i])[0], (*(*_points)[i])[1], (*(*_points)[i])[2]};
		vtkIdType pid = newPoints->InsertNextPoint(coords);
		newVerts->InsertNextCell(1, &pid);
	}

	int nIdx = static_cast<int>(_pnt_idx->size());
	for (int i=0; i<nIdx; i++)
	{
		//size_t idx = (*_pnt_idx)[i];
		//double coords[3] = {(*(*_points)[idx])[0], (*(*_points)[idx])[1], (*(*_points)[idx])[2]};
		//vtkIdType pid = newPoints->InsertNextPoint(coords);
		//newVerts->InsertNextCell(1, &pid);
		vtkIdType id = static_cast<int>((*_pnt_idx)[i]);
		newVerts->InsertNextCell(1, &id);
	}

	if (!_ply_idx->empty())
	{
		for (size_t j=0; j<_ply_idx->size(); j++)
		{
			size_t idx = (*_ply_idx)[j];
			const int nPoints = (*_polylines)[idx]->getNumberOfPoints();
			/*
			for (int i = 0; i < nPoints; i++)
			{
				const GEOLIB::Point* point = (*(*_polylines)[idx])[i];
				const double* coords = point->getData();
				newPoints->InsertNextPoint(coords);
			}
			*/
			// Generate lines
			newLines->InsertNextCell(nPoints);
			//plyIDs->InsertNextValue(j);
			for (int i = 0; i < nPoints; i++)
			{
				int p = (*(*_polylines)[idx]).getPointID(i);
				newLines->InsertCellPoint((*(*_polylines)[idx]).getPointID(i));
			}

		}
	}

	if (!_sfc_idx->empty())
	{
		for (size_t k=0; k<_sfc_idx->size(); k++)
		{

			// benötigte punkte für polys einfugen!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			size_t idx = (*_sfc_idx)[k];
			
			const size_t nTriangles = (*_surfaces)[idx]->getNTriangles();

			for (size_t i = 0; i < nTriangles; i++)
			{
				vtkPolygon* aPolygon = vtkPolygon::New();
				aPolygon->GetPointIds()->SetNumberOfIds(3);

				const GEOLIB::Triangle* triangle = (*(*_surfaces)[idx])[i];
				for (size_t j=0; j<3; j++)
				{
					aPolygon->GetPointIds()->SetId(j, ((*triangle)[j]));
				}
				newPolys->InsertNextCell(aPolygon);
				//sfcIDs->InsertNextValue(count);

				aPolygon->Delete();
			}
		}
	}

	output->SetPoints(newPoints);
	output->SetVerts(newVerts);
	output->SetLines(newLines);
	output->SetPolys(newPolys);

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
