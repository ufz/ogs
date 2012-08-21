/*
 * PetrelInterface.cpp
 *
 *  Created on: Feb 16, 2010
 *      Author: fischeth
 */

#include "PetrelInterface.h"
#include "StringTools.h"
#include <fstream>

namespace FileIO
{
PetrelInterface::PetrelInterface(std::list<std::string> &sfc_fnames,
                                 std::list<std::string> &well_path_fnames,
                                 std::string &unique_model_name, GEOLIB::GEOObjects* geo_obj) :
	_unique_name(unique_model_name), pnt_vec(new std::vector<GEOLIB::Point*>),
	well_vec(new std::vector<GEOLIB::Point*>), ply_vec(new std::vector<
	                                                           GEOLIB::Polyline*>)
{
	for (std::list<std::string>::const_iterator it(sfc_fnames.begin()); it
	     != sfc_fnames.end(); ++it)
	{
		std::cout << "PetrelInterface::PetrelInterface open surface stream from file " <<
		*it
		          << " ... " << std::flush;
		std::ifstream in((*it).c_str());
		if (in)
		{
			std::cout << "done" << std::endl;
			readPetrelSurface(in);
			in.close();
		}
		else
			std::cerr << "error opening stream " << std::endl;
	}

	for (std::list<std::string>::const_iterator it(well_path_fnames.begin()); it
	     != well_path_fnames.end(); ++it)
	{
		std::cout << "PetrelInterface::PetrelInterface open well path stream from file "
		          << *it << " ... " << std::flush;
		std::ifstream in((*it).c_str());
		if (in)
		{
			std::cout << "done" << std::endl;
			readPetrelWellTrace(in);
			in.close();
		}
		else
			std::cerr << "error opening stream " << std::endl;
	}

	// store data in GEOObject
	geo_obj->addPointVec(pnt_vec, _unique_name);
	if (well_vec->size() > 0)
		geo_obj->addStationVec(
		        well_vec,
		        _unique_name);
	if (ply_vec->size () > 0)
		geo_obj->addPolylineVec(ply_vec, _unique_name);
}

PetrelInterface::~PetrelInterface()
{}

void PetrelInterface::readPetrelSurface(std::istream &in)
{
	char buffer[MAX_COLS_PER_ROW];
	in.getline(buffer, MAX_COLS_PER_ROW);
	std::string line(buffer);

	if (line.find("# Petrel Points with attributes") != std::string::npos)
	{
		// read header
		// read Version string
		in.getline(buffer, MAX_COLS_PER_ROW);
		line = buffer;
		// read string BEGIN HEADER
		in.getline(buffer, MAX_COLS_PER_ROW);
		line = buffer;

		in.getline(buffer, MAX_COLS_PER_ROW);
		line = buffer;
		while (line.find("END HEADER") == std::string::npos)
		{
			in.getline(buffer, MAX_COLS_PER_ROW);
			line = buffer;
		}

		// read points
		size_t idx(pnt_vec->size());
		while (in)
		{
			pnt_vec->push_back(new GEOLIB::Point);
			in >> *((*pnt_vec)[idx]);
			if (!in)
			{
				delete (*pnt_vec)[idx];
				pnt_vec->pop_back();
			}
			else
				idx++;
		}
	}
	else
		std::cerr << "error reading petrel points: " << line << std::endl;
}

void PetrelInterface::readPetrelWellTrace(std::istream &in)
{
	char buffer[MAX_COLS_PER_ROW];
	in.getline(buffer, MAX_COLS_PER_ROW);
	std::string line(buffer);

	if (line.find("# WELL TRACE FROM PETREL") != std::string::npos)
	{
		// read header
		// read well name
		in.getline(buffer, MAX_COLS_PER_ROW);
		line = buffer;
		std::list<std::string> str_list(splitString(line, ' '));
		std::list<std::string>::const_iterator it(str_list.begin());
		while (it != str_list.end())
			std::cout << *it++ << " " << std::flush;
		std::cout << std::endl;

		// read well head x coordinate
		in.getline(buffer, MAX_COLS_PER_ROW);
		line = buffer;
		str_list = splitString(line, ' ');
		it = str_list.begin();
		while (it != str_list.end())
			std::cout << *it++ << " " << std::flush;
		std::cout << std::endl;
		it = (str_list.end())--;
		it--;
		char* buf;
		double well_head_x(strtod((*it).c_str(), &buf));

		// read well head y coordinate
		in.getline(buffer, MAX_COLS_PER_ROW);
		line = buffer;
		str_list = splitString(line, ' ');
		it = str_list.begin();
		while (it != str_list.end())
			std::cout << *it++ << " " << std::flush;
		std::cout << std::endl;
		it = (str_list.end())--;
		it--;
		double well_head_y(strtod((*it).c_str(), &buf));

		// read well KB
		in.getline(buffer, MAX_COLS_PER_ROW);
		line = buffer;
		str_list = splitString(line, ' ');
		it = str_list.begin();
		while (it != str_list.end())
			std::cout << *it++ << " " << std::flush;
		std::cout << std::endl;
		it = (str_list.end())--;
		it--;
		double well_kb(strtod((*it).c_str(), &buf));

		std::cout << "PetrelInterface::readPetrelWellTrace: " << well_head_x << "," <<
		well_head_y << "," << well_kb << std::endl;
		well_vec->push_back(
		        static_cast<GEOLIB::StationBorehole*> (new GEOLIB::StationBorehole(
		                                                       well_head_x, well_head_y,
		                                                       well_kb)));

		// read well type
		in.getline(buffer, MAX_COLS_PER_ROW);
		line = buffer;
		str_list = splitString(line, ' ');
		it = str_list.begin();
		while (it != str_list.end())
			std::cout << *it++ << " " << std::flush;
		std::cout << std::endl;
		std::string type(*((str_list.end())--));

		readPetrelWellTraceData(in);
	}
}

void PetrelInterface::readPetrelWellTraceData(std::istream &in)
{
	char buffer[MAX_COLS_PER_ROW];
	in.getline(buffer, MAX_COLS_PER_ROW);
	std::string line(buffer);
	std::list<std::string> str_list;
	std::list<std::string>::const_iterator it;

	// read yet another header lines
	in.getline(buffer, MAX_COLS_PER_ROW);
	line = buffer;
	while (line.substr(0, 1).compare("#") == 0)
	{
		in.getline(buffer, MAX_COLS_PER_ROW);
		line = buffer;
	}

	// read column information
	str_list = splitString(line, ' ');
	it = str_list.begin();
	while (it != str_list.end())
		std::cout << *it++ << " " << std::flush;
	std::cout << std::endl;

	// read points
	double md, x, y, z, tvd, dx, dy, azim, incl, dls;
	in.getline(buffer, MAX_COLS_PER_ROW);
	line = buffer;
	while (in)
	{
		if (line.size() > 1 && line.substr(0, 1).compare("#") != 0)
		{
			std::stringstream stream(line);
			stream >> md;
			stream >> x >> y >> z;
			//			pnt_vec->push_back (new GEOLIB::Point (x,y,z));
			static_cast<GEOLIB::StationBorehole*> ((*well_vec)[well_vec->size()
			                                                   - 1])->addSoilLayer(
			        x,
			        y,
			        z,
			        "unknown");
			stream >> tvd >> dx >> dy >> azim >> incl >> dls;
		}
		in.getline(buffer, MAX_COLS_PER_ROW);
		line = buffer;
	}
}
} // end namespace FileIO
