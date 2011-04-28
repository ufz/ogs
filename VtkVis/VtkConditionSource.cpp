/**
 * \file VtkConditionSource.cpp
 * 2011/03/02 KR Initial implementation
 */

// ** INCLUDES **
#include "VtkConditionSource.h"
#include "AxisAlignedBoundingBox.h"
#include "FEMCondition.h"

#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkSmartPointer.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolygon.h>

#include <vtkLookupTable.h>

vtkStandardNewMacro(VtkConditionSource);
vtkCxxRevisionMacro(VtkConditionSource, "$Revision$");

VtkConditionSource::VtkConditionSource()
: _points(NULL), _cond_vec(NULL)
{
	this->SetNumberOfInputPorts(0);

	const GEOLIB::Color* c = GEOLIB::getRandomColor();
	GetProperties()->SetColor((*c)[0]/255.0,(*c)[1]/255.0,(*c)[2]/255.0);
}

VtkConditionSource::~VtkConditionSource()
{
}

void VtkConditionSource::setData(const std::vector<GEOLIB::Point*>* points, const std::vector<FEMCondition*>* conds)
{
	_points    = points;
	_cond_vec  = conds;
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


	vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkInformation> outInfo = outputVector->GetInformationObject(0);
	vtkSmartPointer<vtkPolyData> output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

	vtkSmartPointer<vtkDoubleArray> scalars = vtkSmartPointer<vtkDoubleArray>::New();
		scalars->SetNumberOfComponents(1);
		scalars->SetName("Scalars");
	//std::map<size_t, size_t> idx_map;

	vtkSmartPointer<vtkCellArray> newVerts = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkCellArray> newLines = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkCellArray> newPolys = vtkSmartPointer<vtkCellArray>::New();

	if (outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()) > 0)
		return 1;

	
	size_t n_pnts = _points->size();
	double value(0.0);
	if (!_cond_vec->empty())
	{
		const std::vector<double> dv = (*_cond_vec)[0]->getDisValue();
		value = dv[dv.size()-1]; // get an existing value for the distribution so scaling on point data will be correct during rendering process!
	}
	
	for (size_t i=0; i<n_pnts; i++)
	{
		double coords[3] = {(*(*_points)[i])[0], (*(*_points)[i])[1], (*(*_points)[i])[2]};
		newPoints->InsertNextPoint(coords);
		scalars->InsertNextValue(value);
	}
	

	size_t nCond = _cond_vec->size();
	for (size_t n=0; n<nCond; n++)
	{
		FiniteElement::DistributionType type = (*_cond_vec)[n]->getProcessDistributionType();
		const std::vector<double> dis_values = (*_cond_vec)[n]->getDisValue();

		if ((*_cond_vec)[n]->getGeoType() == GEOLIB::POINT)
		{
			size_t nPoints = _points->size();
			const GEOLIB::Point* pnt = static_cast<const GEOLIB::Point*>((*_cond_vec)[n]->getGeoObj());
			int id(-1);
			for (size_t i=0; i<nPoints; i++)
			{
				if ((*_points)[i] == pnt) 
				{
					vtkIdType id = static_cast<int>(i);//(this->getIndex(i, newPoints, scalars, idx_map));
					newVerts->InsertNextCell(1, &id);
					if (type == FiniteElement::CONSTANT)
						scalars->SetValue(id, dis_values[0]);
					break;
				}
			}
			if (id==-1) std::cout << "Error in VtkConditionSource::RequestData() - Point object not found ..." << std::endl;
		}
		else if ((*_cond_vec)[n]->getGeoType() == GEOLIB::POLYLINE)
		{
			const GEOLIB::Polyline* ply = static_cast<const GEOLIB::Polyline*>((*_cond_vec)[n]->getGeoObj());
			const int nPoints = ply->getNumberOfPoints();
			newLines->InsertNextCell(nPoints);
			for (int i = 0; i < nPoints; i++)
			{
				size_t pnt_id = ply->getPointID(i);//this->getIndex(ply->getPointID(i), newPoints, scalars, idx_map);
				newLines->InsertCellPoint(pnt_id);

				if (type == FiniteElement::CONSTANT)
					scalars->SetValue(pnt_id, dis_values[0]);
				else if (type == FiniteElement::LINEAR)
				{
					for (size_t j=0; j<dis_values.size(); j+=2)
						if (static_cast<size_t>(dis_values[j]) == pnt_id)
						//if (this->getIndex(static_cast<size_t>(dis_values[j]), newPoints, scalars, idx_map) == pnt_id)
						{
							scalars->SetValue(pnt_id, dis_values[j+1]);
							break;
						}
				}
			}
		}
		else if ((*_cond_vec)[n]->getGeoType() == GEOLIB::SURFACE)
		{
			const GEOLIB::Surface* sfc = static_cast<const GEOLIB::Surface*>((*_cond_vec)[n]->getGeoObj());
			
			const size_t nTriangles = sfc->getNTriangles();

			for (size_t i=0; i<nTriangles; i++)
			{
				vtkPolygon* aPolygon = vtkPolygon::New();
				aPolygon->GetPointIds()->SetNumberOfIds(3);

				const GEOLIB::Triangle* triangle = (*sfc)[i];
				for (size_t j=0; j<3; j++)
				{
					size_t pnt_id = (*triangle)[j];//this->getIndex((*triangle)[j], newPoints, scalars, idx_map);
					aPolygon->GetPointIds()->SetId(j, pnt_id);

					if (type == FiniteElement::CONSTANT)
						scalars->SetValue(pnt_id, dis_values[0]);
					else if (type == FiniteElement::LINEAR)
					{
						for (size_t k=0; k<dis_values.size(); k+=2)
							if (static_cast<size_t>(dis_values[j]) == pnt_id)
							//if (this->getIndex(static_cast<size_t>(dis_values[j]), newPoints, scalars, idx_map) == pnt_id)
							{
								scalars->SetValue(pnt_id, dis_values[j+1]);
								break;
							}
					}
				}
				newPolys->InsertNextCell(aPolygon);

				aPolygon->Delete();
			}
		}
		// draw a bounding box in case of of the conditions is "domain"
		else if ((*_cond_vec)[n]->getGeoType() == GEOLIB::GEODOMAIN)
		{
			GEOLIB::AABB bounding_box (_points);
			std::vector<GEOLIB::Point> box;
			box.push_back(bounding_box.getMinPoint());
			box.push_back(bounding_box.getMaxPoint());

			vtkIdType nPoints = newPoints->GetNumberOfPoints();
			//size_t pnt_idx = _points->size();

			for (size_t i=0; i<8; i++)
			{
				double coords[3] = {box[i%2][0], box[(i>>1)%2][1], box[i>>2][2]};
				newPoints->InsertNextPoint(coords);
				scalars->InsertNextValue(0.0);
				//idx_map.insert( std::pair<size_t,size_t>(pnt_idx+i, nPoints+i));
			}
			
			for (size_t i=0; i<4; i++)
			{
				vtkIdType a[2] = {nPoints+i, nPoints+i+4};
				vtkIdType b[2] = {nPoints+(i*2), nPoints+(i*2+1)};
				vtkIdType c[2] = {nPoints+(static_cast<int>(i/2)*4+(i%2)), nPoints+(static_cast<int>(i/2)*4+(i%2)+2)};
				newLines->InsertNextCell(2, &a[0]);
				newLines->InsertNextCell(2, &b[0]);
				newLines->InsertNextCell(2, &c[0]);
			}
		}
	}


	output->SetPoints(newPoints);
	output->GetPointData()->AddArray(scalars);
	output->GetPointData()->SetActiveScalars("Scalars");
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

/*
size_t VtkConditionSource::getIndex(size_t idx, vtkSmartPointer<vtkPoints> newPoints, vtkSmartPointer<vtkDoubleArray> scalars, std::map<size_t, size_t> &idx_map)
{
	std::map<size_t,size_t>::iterator mapped_index = idx_map.find(idx);
	if (mapped_index != idx_map.end()) return mapped_index->second;

	double coords[3] = {(*(*_points)[idx])[0], (*(*_points)[idx])[1], (*(*_points)[idx])[2]};
	newPoints->InsertNextPoint(coords);
	scalars->InsertNextValue(0.0);
	size_t new_idx = idx_map.size();
	idx_map.insert( std::pair<size_t,size_t>(idx, new_idx) );
	std::cout << idx << ", " << new_idx << std::endl;
	return new_idx;
}
*/
