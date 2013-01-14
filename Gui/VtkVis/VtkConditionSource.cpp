/**
 * \file
 * \author Karsten Rink
 * \date   2011-03-02
 * \brief  Implementation of the VtkConditionSource class.
 *
 * \copyright
 * Copyright (c)  2013, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
// GeoLib
#include "AABB.h"
#include "Color.h"

#include "FEMCondition.h"
#include "VtkConditionSource.h"

#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolygon.h>
#include <vtkSmartPointer.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkProperty.h>

#include <vtkLookupTable.h>

vtkStandardNewMacro(VtkConditionSource);
vtkCxxRevisionMacro(VtkConditionSource, "$Revision$");

VtkConditionSource::VtkConditionSource()
	: _points(NULL), _cond_vec(NULL)
{
	this->SetNumberOfInputPorts(0);

	const GeoLib::Color* c = GeoLib::getRandomColor();
	GetProperties()->SetColor((*c)[0] / 255.0,(*c)[1] / 255.0,(*c)[2] / 255.0);
}

VtkConditionSource::~VtkConditionSource()
{
}

void VtkConditionSource::setData(const std::vector<GeoLib::Point*>* points,
                                 const std::vector<FEMCondition*>* conds)
{
	_removable = false; // From VtkAlgorithmProperties
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

int VtkConditionSource::RequestData( vtkInformation* request,
                                     vtkInformationVector** inputVector,
                                     vtkInformationVector* outputVector )
{
	(void)request;
	(void)inputVector;

	if (this->_points->empty() || this->_cond_vec->empty())
		return 0;

	vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkInformation> outInfo = outputVector->GetInformationObject(0);
	vtkSmartPointer<vtkPolyData> output =
	        vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

	vtkSmartPointer<vtkIdTypeArray> distypes = vtkSmartPointer<vtkIdTypeArray>::New();
	distypes->SetNumberOfComponents(1);
	distypes->SetName("DisTypes");

	vtkSmartPointer<vtkDoubleArray> scalars = vtkSmartPointer<vtkDoubleArray>::New();
	scalars->SetNumberOfComponents(1);
	scalars->SetName("Scalars");
	//std::map<size_t, size_t> idx_map;

	vtkSmartPointer<vtkCellArray> newVerts = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkCellArray> newLines = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkCellArray> newPolys = vtkSmartPointer<vtkCellArray>::New();

	if (outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()) > 0)
		return 1;
	/*
	size_t n_pnts = _points->size();
	double value(-9999);
	if (!_cond_vec->empty())
	{
		const std::vector<double> dv = (*_cond_vec)[0]->getDisValue();
		value = dv[dv.size()-1]; // get an existing value for the distribution so scaling on point data will be correct during rendering process!
	}

	for (size_t i = 0; i < n_pnts; i++)
	{
		double coords[3] = {(*(*_points)[i])[0], (*(*_points)[i])[1], (*(*_points)[i])[2]};
		newPoints->InsertNextPoint(coords);
		distypes->InsertNextValue(0);
		scalars->InsertNextValue(value);
	}
	*/
	vtkIdType pnt_id = 0;
	size_t nCond = _cond_vec->size();
	for (size_t n = 0; n < nCond; n++)
	{
		FiniteElement::DistributionType type = (*_cond_vec)[n]->getProcessDistributionType();
		const std::vector<size_t> dis_nodes = (*_cond_vec)[n]->getDisNodes();
		const std::vector<double> dis_values = (*_cond_vec)[n]->getDisValues();

		vtkIdType dis_type_value(0);
		std::map<FiniteElement::DistributionType, vtkIdType>::const_iterator it(_dis_type_map.find(type));
		if (it == _dis_type_map.end())
		{
			dis_type_value = static_cast<vtkIdType>(_dis_type_map.size());
			_dis_type_map.insert(std::pair<FiniteElement::DistributionType, size_t>(type, dis_type_value));
		}
		else dis_type_value = it->second;

		if ((*_cond_vec)[n]->getGeomType() == GeoLib::POINT)
		{
			/*
			size_t nPoints = _points->size();
			const GeoLib::Point* pnt =
			        static_cast<const GeoLib::Point*>((*_cond_vec)[n]->getGeoObj());
			int id(-1);
			for (size_t i = 0; i < nPoints; i++)
				if ((*_points)[i] == pnt)
				{
					id = static_cast<int>(i); //(this->getIndex(i, newPoints, scalars, idx_map));
					vtkIdType vtk_id = static_cast<vtkIdType>(id);
					*/
					const GeoLib::Point* pnt = static_cast<const GeoLib::Point*>((*_cond_vec)[n]->getGeoObj());
					newPoints->InsertNextPoint(pnt->getCoords());

					newVerts->InsertNextCell(1, &pnt_id);
					if (type == FiniteElement::CONSTANT || type == FiniteElement::CONSTANT_NEUMANN)
						scalars->InsertNextValue(dis_values[0]);
					else scalars->InsertNextValue(0);
					distypes->InsertNextValue(dis_type_value);
					pnt_id++;
			/*
					break;
				}
			if (id == -1)
				std::cout <<
				"Error in VtkConditionSource::RequestData() - Point object not found ..."
				          << std::endl;
			*/
		}
		else if ((*_cond_vec)[n]->getGeomType() == GeoLib::POLYLINE)
		{
			const GeoLib::Polyline* ply = static_cast<const GeoLib::Polyline*>((*_cond_vec)[n]->getGeoObj());
			const int nPoints = ply->getNumberOfPoints();
			newLines->InsertNextCell(nPoints);
			double value (0);
			for (int i = 0; i < nPoints; i++)
			{
				size_t point_index = ply->getPointID(i);

				newPoints->InsertNextPoint((*_points)[point_index]->getCoords());
				newLines->InsertCellPoint(pnt_id);
				distypes->InsertNextValue(dis_type_value);

				if (type == FiniteElement::CONSTANT || type == FiniteElement::CONSTANT_NEUMANN)
					scalars->InsertNextValue(dis_values[0]);
				else if (type == FiniteElement::LINEAR || type == FiniteElement::LINEAR_NEUMANN)
				{
					for (size_t j = 0; j < dis_values.size(); j ++)
					{
						if (static_cast<int>(dis_nodes[j]) == i)
							value = dis_values[j];
					}
					scalars->InsertNextValue(value);
				}
				else
					scalars->InsertNextValue(0);
				pnt_id++;
			}
		}
		else if ((*_cond_vec)[n]->getGeomType() == GeoLib::SURFACE)
		{
			std::vector<int> point_idx_map(_points->size(), -1);

			const GeoLib::Surface* sfc =
			        static_cast<const GeoLib::Surface*>((*_cond_vec)[n]->getGeoObj());

			const size_t nTriangles = sfc->getNTriangles();

			for (size_t i = 0; i < nTriangles; i++)
			{
				vtkPolygon* aPolygon = vtkPolygon::New();
				aPolygon->GetPointIds()->SetNumberOfIds(3);

				const GeoLib::Triangle* triangle = (*sfc)[i];
				for (size_t j = 0; j < 3; j++)
				{
					size_t point_index ((*triangle)[j]);

					if (point_idx_map[point_index] == -1)
					{
						point_idx_map[point_index] = pnt_id;
						newPoints->InsertNextPoint((*_points)[point_index]->getCoords());
						aPolygon->GetPointIds()->SetId(j, pnt_id);
						distypes->InsertNextValue(dis_type_value);

						if (type == FiniteElement::CONSTANT || type == FiniteElement::CONSTANT_NEUMANN)
							scalars->InsertNextValue(dis_values[0]);
						else if (type == FiniteElement::LINEAR || type == FiniteElement::LINEAR_NEUMANN)
						{
							for (size_t k = 0; k < dis_values.size(); k++)
								if (static_cast<size_t>(dis_nodes[j]) == point_index)
								{
									scalars->InsertNextValue(dis_values[j]);
									break;
								}
						}
						else scalars->InsertNextValue(0);
						pnt_id++;
					}
					else
						aPolygon->GetPointIds()->SetId(j, static_cast<vtkIdType>(point_idx_map[point_index]));
				}
				newPolys->InsertNextCell(aPolygon);

				aPolygon->Delete();
			}
		}
		// HACK: this is currently used when applying DIRECT conditions
		else if ((*_cond_vec)[n]->getGeomType() == GeoLib::INVALID)
		{
			size_t nValues = dis_values.size();
			for (size_t i=0; i<nValues; i++)
			{
				//vtkIdType pid = newPoints->InsertNextPoint((*_points)[dis_nodes[i]]->getData());
				vtkIdType pid = newPoints->InsertNextPoint((*_points)[i]->getCoords());
				newVerts->InsertNextCell(1, &pid);
				scalars->InsertNextValue(dis_values[i]);
				distypes->InsertNextValue(dis_type_value);
				pnt_id++;
			}
		}
		// draw a bounding box in case of of the conditions is "domain"
		else if ((*_cond_vec)[n]->getGeomType() == GeoLib::GEODOMAIN)
		{
			GeoLib::AABB<GeoLib::Point> bounding_box (_points->begin(), _points->end());
			std::vector<GeoLib::Point> box;
			box.push_back(bounding_box.getMinPoint());
			box.push_back(bounding_box.getMaxPoint());

			vtkIdType nPoints = newPoints->GetNumberOfPoints();
			//size_t pnt_idx = _points->size();

			for (size_t i = 0; i < 8; i++)
			{
				double coords[3] =
				{box[i % 2][0], box[(i >> 1) % 2][1], box[i >> 2][2]};
				newPoints->InsertNextPoint(coords);
				distypes->InsertNextValue(dis_type_value);
				scalars->InsertNextValue(0.0);
				//idx_map.insert( std::pair<size_t,size_t>(pnt_idx+i, nPoints+i));
			}

			for (vtkIdType i = 0; i < 4; i++)
			{
				vtkIdType a[2] = {nPoints + i, nPoints + i + 4};
				vtkIdType b[2] = {nPoints + (i * 2), nPoints + (i * 2 + 1)};
				vtkIdType c[2] = {nPoints + (static_cast<int>(i / 2) * 4 + (i % 2)), nPoints + (static_cast<int>(i / 2) * 4 + (i % 2) + 2)};
				newLines->InsertNextCell(2, &a[0]);
				newLines->InsertNextCell(2, &b[0]);
				newLines->InsertNextCell(2, &c[0]);
			}
		}
	}

	output->SetPoints(newPoints);
	output->GetPointData()->AddArray(distypes);
	output->GetPointData()->AddArray(scalars);
	output->GetPointData()->SetActiveScalars("Scalars");
	output->SetVerts(newVerts);
	output->SetLines(newLines);
	output->SetPolys(newPolys);

	return 1;
}

int VtkConditionSource::RequestInformation( vtkInformation* request,
                                            vtkInformationVector** inputVector,
                                            vtkInformationVector* outputVector )
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
