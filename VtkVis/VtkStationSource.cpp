/**
 * \file VtkStationSource.cpp
 * 24/02/2010 KR Initial implementation
 *
 */

#include "Station.h"
#include "Color.h"

// ** VTK INCLUDES **
#include "VtkStationSource.h"

#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include "vtkObjectFactory.h"
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkLine.h>
#include <vtkPolyData.h>
#include <vtkCellData.h>

vtkStandardNewMacro(VtkStationSource);
vtkCxxRevisionMacro(VtkStationSource, "$Revision$");

VtkStationSource::VtkStationSource()
: _stations(NULL)
{
	this->SetNumberOfInputPorts(0);

	const GEOLIB::Color* c = GEOLIB::getRandomColor();
	GetProperties()->SetColor((*c)[0]/255.0,(*c)[1]/255.0,(*c)[2]/255.0);
	delete c;
}

VtkStationSource::~VtkStationSource()
{
	std::map<std::string, GEOLIB::Color*>::iterator it;
	for (it = _colorLookupTable.begin(); it != _colorLookupTable.end(); ++it) {
		delete it->second;
	}
}

void VtkStationSource::PrintSelf( ostream& os, vtkIndent indent )
{
	this->Superclass::PrintSelf(os,indent);

	if (_stations->size() == 0)
		return;

	os << indent << "== VtkStationSource ==" << "\n";

	int i = 0;
	for (std::vector<GEOLIB::Point*>::const_iterator it = _stations->begin();
		it != _stations->end(); ++it)
	{
		const double* coords = (*it)->getData();
		os << indent << "Station " << i <<" (" << coords[0] << ", " << coords[1] << ", " << coords[2] << ")\n";
		i++;
	}
}


/// Create 3d Station objects
int VtkStationSource::RequestData( vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector )
{
	(void)request;
	(void)inputVector;

	if (!_stations)
		return 0;
	int nStations = _stations->size();
	if (nStations == 0)
		return 0;

	bool isBorehole = (static_cast<GEOLIB::Station*>((*_stations)[0])->type() == GEOLIB::Station::BOREHOLE) ? true : false;

	vtkSmartPointer<vtkInformation> outInfo = outputVector->GetInformationObject(0);
	vtkSmartPointer<vtkPolyData> output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

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

	vtkSmartPointer<vtkIntArray> strat_ids = vtkSmartPointer<vtkIntArray>::New();
		strat_ids->SetNumberOfComponents(1);
		strat_ids->SetName("Stratigraphies");

	size_t lastMaxIndex(0);
	size_t site_count(0);

	// Generate graphic objects
	for (std::vector<GEOLIB::Point*>::const_iterator it = _stations->begin();
		it != _stations->end(); ++it)
	{
		double coords[3] = { (*(*it))[0], (*(*it))[1], (*(*it))[2] };
		vtkIdType sid = newStations->InsertNextPoint(coords);
		station_ids->InsertNextValue(site_count);

		if (!isBorehole)
		{
			newVerts->InsertNextCell(1, &sid);
		}
		else
		{
			std::vector<GEOLIB::Point*> profile = static_cast<GEOLIB::StationBorehole*>(*it)->getProfile();
			std::vector<std::string> soilNames = static_cast<GEOLIB::StationBorehole*>(*it)->getSoilNames();
			const size_t nLayers = profile.size();

			for (size_t i=1; i<nLayers; i++)
			{
				double* pCoords = const_cast<double*>(profile[i]->getData());
				double loc[3] = { pCoords[0], pCoords[1], pCoords[2] };
				newStations->InsertNextPoint(loc);
				station_ids->InsertNextValue(site_count);

				newLines->InsertNextCell(2);
				newLines->InsertCellPoint(lastMaxIndex);	// start of borehole-layer
				newLines->InsertCellPoint(lastMaxIndex+1);	//end of boreholelayer
				lastMaxIndex++;
				strat_ids->InsertNextValue(this->GetIndexByName(soilNames[i]));
			}
			lastMaxIndex++;
		}
		site_count++;
	}

	output->SetPoints(newStations);

	if (!isBorehole)
	{
		output->SetVerts(newVerts);
		output->GetPointData()->AddArray(station_ids);	
	}
	else
	{
		output->SetLines(newLines);
		output->GetCellData()->AddArray(strat_ids);
		output->GetCellData()->SetActiveAttribute("Stratigraphies", vtkDataSetAttributes::SCALARS);
	}

	return 1;
}

int VtkStationSource::RequestInformation( vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector )
{
	(void)request;
	(void)inputVector;

	vtkInformation* outInfo = outputVector->GetInformationObject(0);
		outInfo->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);

	return 1;
}

void VtkStationSource::SetUserProperty( QString name, QVariant value )
{
	Q_UNUSED(name);
	Q_UNUSED(value);
}

size_t VtkStationSource::GetIndexByName( std::string name )
{
	vtkIdType max_key(0);
	for (std::map<std::string, vtkIdType>::const_iterator it=_id_map.begin(); it != _id_map.end(); ++it)
	{
		if (name.compare(it->first) == 0)
			return it->second;
		if (it->second > max_key) max_key = it->second;
	}
	vtkIdType new_index(max_key+1);
	std::cout << "Key \"" << name << "\" not found in color lookup table..." << std::endl;
	_id_map.insert(std::pair<std::string, vtkIdType>(name, new_index));
	return new_index;
}
