/**
 * \file
 * \author Karsten Rink
 * \date   2010-02-24
 * \brief  Implementation of the VtkStationSource class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** VTK INCLUDES **
#include "VtkStationSource.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "StationBorehole.h"

#include "vtkObjectFactory.h"
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkLine.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkProperty.h>

vtkStandardNewMacro(VtkStationSource);

VtkStationSource::VtkStationSource()
	: _stations(NULL)
{
	_removable = false; // From VtkAlgorithmProperties
	this->SetNumberOfInputPorts(0);

	const GeoLib::Color* c = GeoLib::getRandomColor();
	GetProperties()->SetColor((*c)[0] / 255.0,(*c)[1] / 255.0,(*c)[2] / 255.0);
	delete c;
}

VtkStationSource::~VtkStationSource()
{
	std::map<std::string, GeoLib::Color*>::iterator it;
	for (it = _colorLookupTable.begin(); it != _colorLookupTable.end(); ++it)
		delete it->second;
}

void VtkStationSource::PrintSelf( ostream& os, vtkIndent indent )
{
	this->Superclass::PrintSelf(os,indent);

	if (_stations->size() == 0)
		return;

	os << indent << "== VtkStationSource ==" << "\n";

	int i = 0;
	for (std::vector<GeoLib::Point*>::const_iterator it = _stations->begin();
	     it != _stations->end(); ++it)
	{
		const double* coords = (*it)->getCoords();
		os << indent << "Station " << i << " (" << coords[0] << ", " << coords[1] <<
		", " << coords[2] << ")\n";
		i++;
	}
}

/// Create 3d Station objects
int VtkStationSource::RequestData( vtkInformation* request,
                                   vtkInformationVector** inputVector,
                                   vtkInformationVector* outputVector )
{
	(void)request;
	(void)inputVector;

	if (!_stations)
		return 0;
	size_t nStations = _stations->size();
	if (nStations == 0)
		return 0;

	bool useStationValues(false);
	double sValue=static_cast<GeoLib::Station*>((*_stations)[0])->getStationValue();
	for (size_t i = 1; i < nStations; i++)
		if (static_cast<GeoLib::Station*>((*_stations)[i])->getStationValue() != sValue)
		{
			useStationValues = true;
			break;
		}

	bool isBorehole =
	        (static_cast<GeoLib::Station*>((*_stations)[0])->type() ==
	GeoLib::Station::StationType::BOREHOLE) ? true : false;

	vtkSmartPointer<vtkInformation> outInfo = outputVector->GetInformationObject(0);
	vtkSmartPointer<vtkPolyData> output =
	        vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

	vtkSmartPointer<vtkPoints> newStations = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> newVerts = vtkSmartPointer<vtkCellArray>::New();
	newVerts->Allocate(nStations);

	vtkSmartPointer<vtkCellArray> newLines;

	if (isBorehole)
		newLines = vtkSmartPointer<vtkCellArray>::New();

	if (outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()) > 0)
		return 1;

	vtkSmartPointer<vtkIntArray> station_ids = vtkSmartPointer<vtkIntArray>::New();
	station_ids->SetNumberOfComponents(1);
	station_ids->SetName("SiteIDs");

	vtkSmartPointer<vtkDoubleArray> station_values = vtkSmartPointer<vtkDoubleArray>::New();
	station_values->SetNumberOfComponents(1);
	station_values->SetName("StationValue");

	vtkSmartPointer<vtkIntArray> strat_ids = vtkSmartPointer<vtkIntArray>::New();
	strat_ids->SetNumberOfComponents(1);
	strat_ids->SetName("Stratigraphies");

	size_t lastMaxIndex(0);
	size_t site_count(0);

	// Generate graphic objects
	for (std::vector<GeoLib::Point*>::const_iterator it = _stations->begin();
	     it != _stations->end(); ++it)
	{
		double coords[3] = { (*(*it))[0], (*(*it))[1], (*(*it))[2] };
		vtkIdType sid = newStations->InsertNextPoint(coords);
		station_ids->InsertNextValue(site_count);
		if (useStationValues)
			station_values->InsertNextValue(static_cast<GeoLib::Station*>(*it)->getStationValue());

		if (!isBorehole)
			newVerts->InsertNextCell(1, &sid);
		else
		{
			std::vector<GeoLib::Point*> profile =
			        static_cast<GeoLib::StationBorehole*>(*it)->getProfile();
			std::vector<std::string> soilNames =
			        static_cast<GeoLib::StationBorehole*>(*it)->getSoilNames();
			const size_t nLayers = profile.size();

			for (size_t i = 1; i < nLayers; i++)
			{
				double* pCoords = const_cast<double*>(profile[i]->getCoords());
				double loc[3] = { pCoords[0], pCoords[1], pCoords[2] };
				newStations->InsertNextPoint(loc);
				station_ids->InsertNextValue(site_count);
				newLines->InsertNextCell(2);
				newLines->InsertCellPoint(lastMaxIndex); // start of borehole-layer
				newLines->InsertCellPoint(++lastMaxIndex); //end of boreholelayer
				strat_ids->InsertNextValue(this->GetIndexByName(soilNames[i]));
				if (useStationValues)
					station_values->InsertNextValue(static_cast<GeoLib::Station*>(*it)->getStationValue());
			}
			lastMaxIndex++;
		}
		site_count++;
	}

	output->SetPoints(newStations);

	if (!isBorehole)
	{
		output->SetVerts(newVerts);
		output->GetCellData()->AddArray(station_ids);
		output->GetCellData()->SetActiveAttribute("SiteIDs", vtkDataSetAttributes::SCALARS);
	}
	else
	{
		output->SetLines(newLines);
		//output->GetCellData()->AddArray(station_ids);
		output->GetCellData()->AddArray(strat_ids);
		output->GetCellData()->SetActiveAttribute("Stratigraphies", vtkDataSetAttributes::SCALARS);
	}
	if (useStationValues)
		output->GetPointData()->AddArray(station_values);

	output->Squeeze();

	return 1;
}

int VtkStationSource::RequestInformation( vtkInformation* /*request*/,
                                          vtkInformationVector** /*inputVector*/,
                                          vtkInformationVector* /*outputVector*/ )
{
	return 1;
}

void VtkStationSource::SetUserProperty( QString name, QVariant value )
{
	Q_UNUSED(name);
	Q_UNUSED(value);
}

size_t VtkStationSource::GetIndexByName( std::string const& name )
{
	vtkIdType max_key(0);
	for (std::map<std::string, vtkIdType>::const_iterator it = _id_map.begin();
	     it != _id_map.end(); ++it)
	{
		if (name.compare(it->first) == 0)
			return it->second;
		if (it->second > max_key)
			max_key = it->second;
	}

	vtkIdType new_index = (_id_map.empty()) ? 0 : (max_key+1);
	INFO("Key \"%s\" has been assigned index %d.", name.c_str(), new_index);
	_id_map.insert(std::pair<std::string, vtkIdType>(name, new_index));
	return new_index;
}

