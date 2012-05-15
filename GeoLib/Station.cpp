/**
 * \file Station.cpp
 * KR Initial implementation
 */

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
// Base
#include "DateTools.h"
#include "StringTools.h"
// GeoLib
#include "Station.h"

namespace GeoLib
{
Station::Station(double x, double y, double z, std::string name) :
	Point (x,y,z), _name(name), _type(Station::STATION), _station_value(0.0)
{
	addProperty("x", &getX, &Station::setX);
	addProperty("y", &getY, &Station::setY);
	addProperty("z", &getZ, &Station::setZ);
}

Station::Station(Point* coords, std::string name) :
	Point (*coords), _name(name), _type(Station::STATION), _station_value(0.0)
{
	addProperty("x", &getX, &Station::setX);
	addProperty("y", &getY, &Station::setY);
	addProperty("z", &getZ, &Station::setZ);
}

Station::Station(Station const& src) :
	Point(src.getCoords()), _name(src._name), _type(src._type),
	_station_value(src._station_value)
{
	addProperty("x", &getX, &Station::setX);
	addProperty("y", &getY, &Station::setY);
	addProperty("z", &getZ, &Station::setZ);
}

void Station::addProperty(std::string pname, double (* getFct)(void*), void (* set)(void*, double))
{
	STNProperty p;
	p.name = pname;
	p.get = getFct;
	p.set = set;
	_properties.push_back(p);
}

Station::~Station()
{
}

Station* Station::createStation(const std::string & line)
{
	std::list<std::string>::const_iterator it;
	Station* station = new Station();
	std::list<std::string> fields = splitString(line, '\t');

	if (fields.size() >= 3)
	{
		it = fields.begin();
		station->_name  = *it;
		(*station)[0]     = strtod((replaceString(",", ".", *(++it))).c_str(), NULL);
		(*station)[1]     = strtod((replaceString(",", ".", *(++it))).c_str(), NULL);
		if (++it != fields.end())
			(*station)[2] = strtod((replaceString(",", ".", *it)).c_str(), NULL);
	}
	else
	{
		std::cout << "Station::createStation() - Unexpected file format..." << std::endl;
		delete station;
		return NULL;
	}
	return station;
}

Station* Station::createStation(const std::string &name, double x, double y, double z)
{
	Station* station = new Station();
	station->_name = name;
	(*station)[0] = x;
	(*station)[1] = y;
	(*station)[2] = z;
	return station;
}

const std::map<std::string, double> Station::getProperties()
{
	std::map<std::string, double> propertyMap;

	for (int i = 0; i < static_cast<int>(_properties.size()); i++)
	{
		double (* getFct)(void*) = _properties[i].get;
		//setFct set = _properties[i].set;
		propertyMap[_properties[i].name] = (*getFct)((void*)this);
	}

	return propertyMap;
}

bool Station::inSelection(const std::vector<PropertyBounds> &bounds)
{
	double value;
	for (size_t i = 0; i < bounds.size(); i++)
	{
		for (size_t j = 0; j < _properties.size(); j++)
			if (_properties[j].name.compare(bounds[i].getName()) == 0)
			{
				double (* get)(void*) = _properties[j].get;
				value = (*get)((void*)this);
				if (!(value >= bounds[i].getMin() && value <= bounds[i].getMax()))
					return false;
			}
	}
	return true;
}

////////////////////////
// The Borehole class //
////////////////////////

StationBorehole::StationBorehole(double x, double y, double z) :
	Station (x,y,z), _zCoord(0), _depth(0), _date(0)
{
	_type = Station::BOREHOLE;
	addProperty("date", &StationBorehole::getDate, &StationBorehole::setDate);
	addProperty("depth", &StationBorehole::getDepth, &StationBorehole::setDepth);

	// add first point of borehole
	_profilePntVec.push_back(this);
	_soilName.push_back("");
}

StationBorehole::~StationBorehole(void)
{
	// deletes profile vector of borehole, starting at layer 1
	// the first point is NOT deleted as it points to the station object itself
	for (size_t k(1); k < _profilePntVec.size(); k++)
		delete _profilePntVec[k];
}

int StationBorehole::find(const std::string &str)
{
	size_t size = _soilName.size();
	for (size_t i = 0; i < size; i++)
		if (_soilName[i].find(str) == 0)
			return 1;
	return 0;
}

int StationBorehole::readStratigraphyFile(const std::string &path,
                                          std::vector<std::list<std::string> > &data)
{
	std::string line;
	std::ifstream in( path.c_str() );

	if (!in.is_open())
	{
		std::cout << "StationBorehole::readStratigraphyFile() - Could not open file..." <<
		std::endl;
		return 0;
	}

	while ( getline(in, line) )
	{
		std::list<std::string> fields = splitString(line, '\t');
		data.push_back(fields);
	}

	in.close();

	return 1;
}

int StationBorehole::addStratigraphy(const std::string &path, StationBorehole* borehole)
{
	std::vector<std::list<std::string> > data;
	if (readStratigraphyFile(path, data))
	{
		size_t size = data.size();
		for (size_t i = 0; i < size; i++)
			addLayer(data[i], borehole);

		// check if a layer is missing
		size = borehole->_soilName.size();
		std::cout << "StationBorehole::addStratigraphy ToDo" << std::endl;
		//	for (size_t i=0; i<size; i++)
		//	{
		//		if ((borehole->_soilLayerThickness[i] == -1) ||(borehole->_soilName[i].compare("") == 0))
		//		{
		//			borehole->_soilLayerThickness.clear();
		//			borehole->_soilName.clear();
		//
		//			cout << "StationBorehole::addStratigraphy() - Profile incomplete (Borehole " << borehole->_name << ", Layer " << (i+1) << " missing).\n";
		//
		//			return 0;
		//		}
		//	}
	}
	else
		borehole->addSoilLayer(borehole->getDepth(), "depth");

	return 1;
}

int StationBorehole::addLayer(std::list<std::string> fields, StationBorehole* borehole)
{
	if (fields.size() >= 4) /* check if there are enough fields to create a borehole object */
	{
		if (fields.front().compare(borehole->_name) == 0) /* check if the name of the borehole matches the name in the data */
		{
			fields.pop_front();

			// int layer = atoi(fields.front().c_str());
			fields.pop_front();

			std::cerr << "StationBorehole::addLayer - assuming correct order" <<
			std::endl;
			double thickness(strtod(replaceString(",", ".", fields.front()).c_str(), 0));
			fields.pop_front();
			borehole->addSoilLayer(thickness, fields.front());
		}
	}
	else
	{
		std::cout
		<< "StationBorehole::addLayer() - Unexpected file format (Borehole "
		<< borehole->_name << ")..." << std::endl;
		return 0;
	}
	return 1;
}

int StationBorehole::addStratigraphy(const std::vector<GeoLib::Point*> &profile, const std::vector<std::string> soil_names)
{
	if (((profile.size()-1) == soil_names.size()) && (soil_names.size()>0))
	{
		this->_profilePntVec.push_back(profile[0]);
		size_t nLayers = soil_names.size();
		for (size_t i=0; i<nLayers; i++)
		{
			this->_profilePntVec.push_back(profile[i+1]);
			this->_soilName.push_back(soil_names[i]);
		}
		return 1;
	}
	
	std::cout << "Error in StationBorehole::addStratigraphy() - Length of parameter vectors does not match." << std::endl;
	return 0;
}

int StationBorehole::addStratigraphies(const std::string &path, std::vector<Point*>* boreholes)
{
	std::vector<std::list<std::string> > data;

	if (readStratigraphyFile(path, data))
	{
		std::string name;

		size_t it = 0;
		size_t nBoreholes = data.size();
		for (size_t i = 0; i < nBoreholes; i++)
		{
			std::list<std::string> fields = data[i];

			if (fields.size() >= 4)
			{
				name = static_cast<StationBorehole*>((*boreholes)[it])->_name;
				if ( fields.front().compare(name) != 0 )
					if (it < boreholes->size() - 1)
						it++;

				fields.pop_front();
				//the method just assumes that layers are read in correct order
				fields.pop_front();
				double thickness (strtod(replaceString(",", ".",
				                                       fields.front()).c_str(), 0));
				fields.pop_front();
				std::string soil_name (fields.front());
				fields.pop_front();
				static_cast<StationBorehole*>((*boreholes)[it])->addSoilLayer(
				        thickness,
				        soil_name);
			}
			else
				std::cout <<
				"StationBorehole::addStratigraphies() - Unexpected file format..."
				          << std::endl;
				//return 0;
		}
	}
	else
		createSurrogateStratigraphies(boreholes);

	return 1;
}

StationBorehole* StationBorehole::createStation(const std::string &line)
{
	StationBorehole* borehole = new StationBorehole();
	std::list<std::string> fields = splitString(line, '\t');

	if (fields.size()      >= 5)
	{
		borehole->_name     = fields.front();
		fields.pop_front();
		(*borehole)[0]      = strtod((replaceString(",", ".", fields.front())).c_str(), NULL);
		fields.pop_front();
		(*borehole)[1]      = strtod((replaceString(",", ".", fields.front())).c_str(), NULL);
		fields.pop_front();
		(*borehole)[2]      = strtod((replaceString(",", ".", fields.front())).c_str(), NULL);
		fields.pop_front();
		borehole->_depth    = strtod((replaceString(",", ".", fields.front())).c_str(), NULL);
		fields.pop_front();
		if (fields.empty())
			borehole->_date = 0;
		else
		{
			borehole->_date = BaseLib::strDate2double(fields.front());
			fields.pop_front();
		}
	}
	else
	{
		std::cout << "Station::createStation() - Unexpected file format..." << std::endl;
		delete borehole;
		return NULL;
	}
	return borehole;
}

StationBorehole* StationBorehole::createStation(const std::string &name,
                                                double x,
                                                double y,
                                                double z,
                                                double depth,
                                                std::string date)
{
	StationBorehole* station = new StationBorehole();
	station->_name  = name;
	(*station)[0]   = x;
	(*station)[1]   = y;
	(*station)[2]   = z;
	station->_depth = depth;
	if (date.compare("0000-00-00") != 0)
		station->_date  = BaseLib::xmlDate2double(date);
	return station;
}

void StationBorehole::createSurrogateStratigraphies(std::vector<Point*>* boreholes)
{
	size_t nBoreholes = boreholes->size();
	for (size_t i = 0; i < nBoreholes; i++)
	{
		StationBorehole* bore = static_cast<StationBorehole*>((*boreholes)[i]);
		bore->addSoilLayer(bore->getDepth(), "depth");
	}
}

void StationBorehole::addSoilLayer ( double thickness, const std::string &soil_name)
{
	/*
	   // TF - Altmark
	   if (_profilePntVec.empty())
	    addSoilLayer ((*this)[0], (*this)[1], (*this)[2]-thickness, soil_name);
	   else {
	    size_t idx (_profilePntVec.size());
	    // read coordinates from last above
	    double x((*_profilePntVec[idx-1])[0]);
	    double y((*_profilePntVec[idx-1])[1]);
	    double z((*_profilePntVec[idx-1])[2]-thickness);
	    addSoilLayer (x, y, z, soil_name);
	   }
	 */

	// KR - Bode
	if (_profilePntVec.empty())
		addSoilLayer ((*this)[0], (*this)[1], (*this)[2], "");

	size_t idx (_profilePntVec.size());
	double x((*_profilePntVec[idx - 1])[0]);
	double y((*_profilePntVec[idx - 1])[1]);
	double z((*_profilePntVec[0])[2] - thickness);
	addSoilLayer (x, y, z, soil_name);
}

void StationBorehole::addSoilLayer ( double x, double y, double z, const std::string &soil_name)
{
	_profilePntVec.push_back (new Point (x, y, z));
	_soilName.push_back(soil_name);
}
} // namespace
