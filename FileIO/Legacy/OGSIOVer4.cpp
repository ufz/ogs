/**
 * \file
 * \author Thomas Fischer
 * \date   2010-01-14
 * \brief  Implementation of the OGSIOVer4 class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <iomanip>
#include <sstream>

// ThirdParty/logog
#include "logog/include/logog.hpp"

// FileIO
#include "MeshIO/GMSHInterface.h"
#include "OGSIOVer4.h"

// BaseLib
#include "FileTools.h"
#include "StringTools.h"
#include "quicksort.h"

// GeoLib
#include "GEOObjects.h"
#include "Point.h"
#include "Polygon.h"
#include "Polyline.h"
#include "SimplePolygonTree.h"
#include "Surface.h"
#include "Triangle.h"

// for tests only
#include "PointVec.h"

// MathLib
#include "AnalyticalGeometry.h"
#include "EarClippingTriangulation.h"

using namespace GeoLib;

namespace FileIO
{
/**************************************************************************
   GeoLib- Funktion: readPoints
   Aufgabe: Lesen der GLI Points und schreiben in einen Vector
   08/2005 CC Implementation
   01/2010 TF big modifications
**************************************************************************/
/** reads the points inclusive their names from input stream in
 * using the OGS-4 file format */
std::string readPoints(std::istream &in, std::vector<Point*>* pnt_vec,
                       bool &zero_based_indexing, std::map<std::string,size_t>* pnt_id_name_map)
{
	std::string line;
	size_t cnt(0);

	getline(in, line);
	// geometric key words start with the hash #
	// while not found a new key word do ...
	while (line.find("#") == std::string::npos && !in.eof() && !in.fail())
	{
		// read id and point coordinates
		std::stringstream inss(line);
		size_t id;
		double x, y, z;
		inss >> id >> x >> y >> z;
		if (!inss.fail ())
		{
			if (cnt == 0)
			{
				if (id == 0)
					zero_based_indexing = true;
				else
					zero_based_indexing = false;
			}
			pnt_vec->push_back(new Point(x, y, z));

			// read mesh density
			if (line.find("$MD") != std::string::npos)
			{
				double mesh_density;
				size_t pos1(line.find_first_of("M"));
				inss.str(line.substr(pos1 + 2, std::string::npos));
				inss >> mesh_density;
			}

			// read name of point
			size_t pos (line.find("$NAME"));
			if (pos != std::string::npos) //OK
			{
				size_t end_pos ((line.substr (pos + 6)).find(" "));
				if (end_pos != std::string::npos)
					(*pnt_id_name_map)[line.substr (pos + 6, end_pos)] = id;
				else
					(*pnt_id_name_map)[line.substr (pos + 6)] = id;
			}

			size_t id_pos (line.find("$ID"));
			if (id_pos != std::string::npos)
				WARN("readPoints(): found tag $ID - please use tag $NAME for reading point names in point %d.", cnt);
			cnt++;
		}
		getline(in, line);
	}

	return line;
}

/** reads points from a vector */
void readPolylinePointVector(const std::string &fname,
                             std::vector<Point*>& pnt_vec,
                             GeoLib::Polyline* ply,
                             const std::string &path,
                             std::vector<std::string> &errors)
{
	// open file
	std::ifstream in((path + fname).c_str());
	if (!in) {
		WARN("readPolylinePointVector(): error opening stream from %s", fname.c_str());
		errors.push_back ("[readPolylinePointVector] error opening stream from " + fname);
		return;
	}

	double x, y, z;
	while (in)
	{
		in >> x >> y >> z;
		size_t pnt_id(pnt_vec.size());
		pnt_vec.push_back(new Point(x, y, z));
		ply->addPoint(pnt_id);
	}
}

/**************************************************************************
   GeoLib-Method: Read
   Task: Read polyline data from file
   Programing:
   03/2004 CC Implementation
   09/2004 OK file path for PLY files
   07/2005 CC PLY id
   08/2005 CC parameter
   09/2005 CC itoa - convert integer to string
   01/2010 TF cleaned method from unused variables
**************************************************************************/
/** read a single Polyline from stream in into the ply_vec-vector */
std::string readPolyline(std::istream &in,
                         std::vector<GeoLib::Polyline*>* ply_vec,
                         std::map<std::string,size_t>& ply_vec_names,
                         std::vector<Point*>& pnt_vec,
                         bool zero_based_indexing,
                         const std::vector<size_t>& pnt_id_map,
                         const std::string &path,
                         std::vector<std::string> &errors)
{
	std::string line, name_of_ply;
	GeoLib::Polyline* ply(new GeoLib::Polyline(pnt_vec));
	size_t type = 2; // need an initial value

	// Schleife ueber alle Phasen bzw. Komponenten
	do {
		in >> line;
		if (line.find("$ID") != std::string::npos) // subkeyword found CC
			in >> line; // read value
			//			id = strtol(line_string.data(), NULL, 0);
		//....................................................................
		if (line.find("$NAME") != std::string::npos) // subkeyword found
		{
			in >> line;
			name_of_ply = line.substr(0); // read value
		}
		//....................................................................
		if (line.find("$TYPE") != std::string::npos) // subkeyword found
		{
			in >> line; // read value
			type = static_cast<size_t> (strtol(line.c_str(), NULL, 0));
		}
		//....................................................................
		if (line.find("$EPSILON") != std::string::npos) // subkeyword found
			in >> line; // read value
		//....................................................................
		if (line.find("$MAT_GROUP") != std::string::npos) // subkeyword found
			in >> line; // read value
		//....................................................................
		if (line.find("$POINTS") != std::string::npos) // subkeyword found
		{ // read the point ids
			in >> line;
			if (type != 100)
				while (!in.eof() && line.size() != 0
				       && (line.find("#") == std::string::npos)
				       && (line.find("$") == std::string::npos))
				{
					size_t pnt_id(BaseLib::str2number<size_t> (line));
					if (!zero_based_indexing)
						pnt_id--;  // one based indexing
					size_t ply_size (ply->getNumberOfPoints());
					if (ply_size > 0)
					{
						if (ply->getPointID (ply_size - 1) != pnt_id_map[pnt_id])
							ply->addPoint(pnt_id_map[pnt_id]);
					}
					else
						ply->addPoint(pnt_id_map[pnt_id]);
					in >> line;
				}
			else {
				WARN("readPolyline(): polyline is an arc *** reading not implemented");
				errors.push_back ("[readPolyline] reading polyline as an arc is not implemented");
			}
			// empty line or the keyword or subkeyword or end of file
		}
		//....................................................................
		if (line.find("$POINT_VECTOR") != std::string::npos) // subkeyword found
		{
			in >> line; // read file name
			line = path + line;
			readPolylinePointVector(line, pnt_vec, ply, path, errors);
		} // subkeyword found
	} while (line.find("#") == std::string::npos && line.size() != 0 && in);

	if (type != 100)
	{
		ply_vec_names.insert (std::pair<std::string,size_t>(name_of_ply, ply_vec->size()));
		ply_vec->push_back(ply);
	}

	return line;
}

/**************************************************************************
   GEOLib-Function:
   Task: polyline read function
   Programming:
   03/2004 CC Implementation
   05/2004 CC Modification
   04/2005 CC Modification calculate the minimal distance between points reference for
   mesh density of line element calculation
   07/2005 CC read ID of polyline
   08/2005 CC parameter
   01/2010 TF changed signature of function
**************************************************************************/
/** reads polylines */
std::string readPolylines(std::istream &in, std::vector<GeoLib::Polyline*>* ply_vec,
                          std::map<std::string,size_t>& ply_vec_names, std::vector<Point*>& pnt_vec,
                          bool zero_based_indexing, const std::vector<size_t>& pnt_id_map,
                          const std::string &path, std::vector<std::string>& errors)
{
	if (!in) {
		WARN("readPolylines(): input stream error.");
		return std::string("");
	}
	std::string tag("#POLYLINE");

	while (!in.eof() && tag.find("#POLYLINE") != std::string::npos)
		tag = readPolyline(in, ply_vec, ply_vec_names, pnt_vec,
		                   zero_based_indexing, pnt_id_map, path, errors);

	return tag;
}

void readTINFile(const std::string &fname, Surface* sfc,
                 std::vector<Point*> &pnt_vec, std::vector<std::string>& errors)
{
	// open file
	std::ifstream in(fname.c_str());
	if (!in) {
		WARN("readTINFile(): could not open stream from %s.", fname.c_str());
		errors.push_back ("readTINFile error opening stream from " + fname);
		return;
	}

	size_t id;
	double x, y, z;
	while (in)
	{
		// read id
		in >> id;
		// determine size
		size_t pnt_pos(pnt_vec.size());
		// read first point
		in >> x >> y >> z;
		pnt_vec.push_back(new Point(x, y, z));
		// read second point
		in >> x >> y >> z;
		pnt_vec.push_back(new Point(x, y, z));
		// read third point
		in >> x >> y >> z;
		pnt_vec.push_back(new Point(x, y, z));
		// create new Triangle
		sfc->addTriangle(pnt_pos, pnt_pos + 1, pnt_pos + 2);
	}
}

/**************************************************************************
   GeoLib-Method: readSurface
   Task: Read surface data from input stream
   Programing:
   03/2004 OK Implementation
   05/2005 OK EPSILON
   09/2005 CC file_path_base
   01/2010 TF signatur modification, reimplementation
**************************************************************************/
/** read a single Surface */
std::string readSurface(std::istream &in,
                        std::vector<GeoLib::Polygon*> &polygon_vec,
                        std::vector<Surface*> &sfc_vec,
                        std::map<std::string,size_t>& sfc_names,
                        const std::vector<GeoLib::Polyline*> &ply_vec,
                        const std::map<std::string, size_t>& ply_vec_names,
                        std::vector<Point*> &pnt_vec,
                        std::string const& path, std::vector<std::string>& errors)
{
	std::string line;
	Surface* sfc(NULL);

	int type (-1);
	std::string name;
	size_t ply_id (0); // std::numeric_limits<size_t>::max());

	do {
		in >> line;
		if (line.find("$ID") != std::string::npos) // subkeyword found CC
			in >> line; // read value
			//			id = strtol(line_string.data(), NULL, 0);
		//....................................................................
		if (line.find("$NAME") != std::string::npos) // subkeyword found
		{
			in >> line; // read value
			name = line.substr(0);
		}
		//....................................................................
		if (line.find("$TYPE") != std::string::npos) // subkeyword found
		{
			in >> line; // read value
			type = strtol(line.c_str(), NULL, 0);
		}
		//....................................................................
		if (line.find("$EPSILON") != std::string::npos) // subkeyword found
			in >> line; // read value
		//....................................................................
		if (line.find("$TIN") != std::string::npos) // subkeyword found
		{
			in >> line; // read value (file name)
			line = path + line;
			sfc = new Surface(pnt_vec);

			readTINFile(line, sfc, pnt_vec, errors);
			if (sfc->getNTriangles() == 0) {
				delete sfc;
				sfc = NULL;
			}
		}
		//....................................................................
		if (line.find("$MAT_GROUP") != std::string::npos) // subkeyword found
			in >> line; // read value
		//....................................................................
		if (line.find("$POLYLINES") != std::string::npos) // subkeyword found
		{ // read the name of the polyline(s)
			in >> line;
			while (!in.eof() && line.size() != 0
			       && (line.find("#") == std::string::npos)
			       && (line.find("$") == std::string::npos))
			{
				// we did read the name of a polyline -> search the id for polyline
				std::map<std::string,size_t>::const_iterator it (ply_vec_names.find (
				                                                         line));
				if (it != ply_vec_names.end())
					ply_id = it->second;
				else
					ply_id = ply_vec.size();

				if (ply_id == ply_vec.size()) {
					WARN("readSurface(): polyline for surface not found!");
					errors.push_back("[readSurface] polyline for surface not found!");
				} else {
					if (type == 3) {
						WARN("readSurface(): surface type 3: flat surface with any normal direction - reading not implemented.");
						errors.push_back("[readSurface] surface type 3: flat surface with any normal direction - reading not implemented");
					}
					if (type == 2) {
						WARN("readSurface(): vertical surface (type 2) - reading not implemented");
						errors.push_back("[readSurface] vertical surface (type 2) - reading not implemented");
					}
				}
				in >> line;
			}
			// empty line or a keyword is found
		}
	} while (line.find("#") == std::string::npos && line.size() != 0 && in);

	if (!name.empty())
		sfc_names.insert(std::pair<std::string,size_t>(name,sfc_vec.size()));

	if (sfc)
		// surface create by TIN
		sfc_vec.push_back (sfc);
	else
    {
        // surface created by polygon
        if (ply_id != std::numeric_limits<size_t>::max() && ply_id != ply_vec.size())
        {
            if (ply_vec[ply_id]->isClosed())
            {
                polygon_vec.push_back (new GeoLib::Polygon (*(ply_vec[ply_id]), true));
            }
            else
            {
            	WARN("readSurface(): cannot create surface %s from polyline %d since polyline is not closed.",  name.c_str(), ply_id);
            }
        }
    }

	return line;
}

/**************************************************************************
   GEOLib-Method:
   Task: Surface read function
   Programming:
   03/2004 OK Implementation
   05/2004 CC Modification
   01/2010 TF changed signature of function, big modifications
**************************************************************************/
std::string readSurfaces(std::istream &in,
                         std::vector<Surface*> &sfc_vec,
                         std::map<std::string, size_t>& sfc_names,
                         const std::vector<GeoLib::Polyline*> &ply_vec,
                         const std::map<std::string,size_t>& ply_vec_names,
                         std::vector<Point*> &pnt_vec,
                         const std::string &path, std::vector<std::string>& errors)
{
	if (!in.good())
	{
		WARN("readSurfaces(): input stream error.");
		return std::string("");
	}
	std::string tag("#SURFACE");

	std::vector<GeoLib::Polygon*> polygon_vec;

	while (!in.eof() && tag.find("#SURFACE") != std::string::npos)
	{
		size_t n_polygons (polygon_vec.size());
		tag = readSurface(in,
		                  polygon_vec,
		                  sfc_vec,
		                  sfc_names,
		                  ply_vec,
		                  ply_vec_names,
		                  pnt_vec,
		                  path,
		                  errors);
		if (n_polygons < polygon_vec.size())
		{
			// subdivide polygon in simple polygons
			GeoLib::Surface* sfc(GeoLib::Surface::createSurface(
			                             *(dynamic_cast<GeoLib::Polyline*> (polygon_vec
			                                                                [
			                                                                        polygon_vec
			                                                                        .
			                                                                        size() - 1]))));
			sfc_vec.push_back(sfc);
		}
	}
	for (size_t k(0); k < polygon_vec.size(); k++)
		delete polygon_vec[k];

	return tag;
}

bool readGLIFileV4(const std::string& fname,
                   GEOObjects* geo,
                   std::string& unique_name,
                   std::vector<std::string>& errors)
{
	INFO("GeoLib::readGLIFile(): open stream from file %s.", fname.c_str());
	std::ifstream in(fname.c_str());
	if (!in) {
		WARN("GeoLib::readGLIFile(): could not open file %s.", fname.c_str());
		errors.push_back("[readGLIFileV4] error opening stream from " + fname);
		return false;
	}
	INFO("GeoLib::readGLIFile(): \t done.");

	std::string tag;
	while (tag.find("#POINTS") == std::string::npos && !in.eof())
		getline (in, tag);

	// read names of points into vector of strings
	std::map<std::string,size_t>* pnt_id_names_map (new std::map<std::string,size_t>);
	bool zero_based_idx(true);
	std::vector<Point*>* pnt_vec(new std::vector<Point*>);
	INFO("GeoLib::readGLIFile(): read points from stream.");
	tag = readPoints(in, pnt_vec, zero_based_idx, pnt_id_names_map);
	INFO("GeoLib::readGLIFile(): \t ok, %d points read.", pnt_vec->size());

	unique_name = BaseLib::extractBaseName(fname);
	if (!pnt_vec->empty())
		geo->addPointVec(pnt_vec, unique_name, pnt_id_names_map);  // KR: insert into GEOObjects if not empty

	// extract path for reading external files
	const std::string path = BaseLib::extractPath(fname);

	// read names of plys into temporary string-vec
	std::map<std::string,size_t>* ply_names (new std::map<std::string,size_t>);
	std::vector<GeoLib::Polyline*>* ply_vec(new std::vector<GeoLib::Polyline*>);
	if (tag.find("#POLYLINE") != std::string::npos && in)
	{
		INFO("GeoLib::readGLIFile(): read polylines from stream.");
		tag = readPolylines(in, ply_vec, *ply_names, *pnt_vec,
		                    zero_based_idx, geo->getPointVecObj(
		                            unique_name)->getIDMap(), path, errors);
		INFO("GeoLib::readGLIFile(): \t ok, %d polylines read.", ply_vec->size());
	}
	else
		INFO("GeoLib::readGLIFile(): tag #POLYLINE not found.");

	std::vector<Surface*>* sfc_vec(new std::vector<Surface*>);
	std::map<std::string,size_t>* sfc_names (new std::map<std::string,size_t>);
	if (tag.find("#SURFACE") != std::string::npos && in)
	{
		INFO("GeoLib::readGLIFile(): read surfaces from stream.");
		tag = readSurfaces(in,
		                   *sfc_vec,
		                   *sfc_names,
		                   *ply_vec,
		                   *ply_names,
		                   *pnt_vec,
		                   path,
		                   errors);
		INFO("GeoLib::readGLIFile(): \tok, %d surfaces read.", sfc_vec->size());
	}
	else
		INFO("GeoLib::readGLIFile(): tag #SURFACE not found.");

	in.close();

	if (!ply_vec->empty())
		geo->addPolylineVec(ply_vec, unique_name, ply_names);  // KR: insert into GEOObjects if not empty
	else
		delete ply_vec;
	if (!sfc_vec->empty())
		geo->addSurfaceVec(sfc_vec, unique_name, sfc_names);  // KR: insert into GEOObjects if not empty
	else
		delete sfc_vec;

	if (errors.empty())
		return true;
	else
		return false;
}

void writeGLIFileV4 (const std::string& fname,
                     const std::string& geo_name,
                     const GeoLib::GEOObjects& geo)
{
	GeoLib::PointVec const* const pnt_vec(geo.getPointVecObj(geo_name));
	std::vector<GeoLib::Point*> const* const pnts (pnt_vec->getVector());
	std::ofstream os (fname.c_str());
	if (pnts) {
		std::string pnt_name;
		const size_t n_pnts(pnts->size());
		INFO("GeoLib::writeGLIFileV4(): writing %d points to file %s.", n_pnts, fname.c_str());
		os << "#POINTS" << std::endl;
		os.precision (20);
		for (size_t k(0); k < n_pnts; k++) {
			os << k << " " << *((*pnts)[k]) << std::flush;
			if (pnt_vec->getNameOfElementByID(k, pnt_name)) {
				os << " $NAME " << pnt_name << std::flush;
			}
			os << std::endl;
		}
	}

	const GeoLib::PolylineVec* plys_vec (geo.getPolylineVecObj (geo_name));
	if (plys_vec)
	{
		const std::vector<GeoLib::Polyline*>* plys (plys_vec->getVector());
		INFO("GeoLib::writeGLIFileV4(): %d polylines to file %s.",
		     plys->size (), fname.c_str());
		for (size_t k(0); k < plys->size(); k++)
		{
			os << "#POLYLINE" << std::endl;
			std::string polyline_name;
			plys_vec->getNameOfElement((*plys)[k], polyline_name);
			os << " $NAME " << std::endl << "  " << polyline_name << std::endl;
			os << " $POINTS" << std::endl;
			for (size_t j(0); j < (*plys)[k]->getNumberOfPoints(); j++)
				os << "  " << ((*plys)[k])->getPointID(j) << std::endl;
		}
	}

	if (plys_vec)
	{
		const std::vector<GeoLib::Polyline*>* plys(plys_vec->getVector());
		INFO("GeoLib::writeGLIFileV4(): write closed polylines as surfaces to file %s.",
		     fname.c_str());
		for (size_t k(0); k < plys->size(); k++)
			if ((*plys)[k]->isClosed())
			{
				os << "#SURFACE" << std::endl;
				os << " $NAME " << std::endl << "  " << k << std::endl; //plys_vec->getNameOfElement ((*plys)[k]) << std::endl;
				os << " $TYPE " << std::endl << "  0" << std::endl;
				os << " $POLYLINES" << std::endl << "  " << k << std::endl; //plys_vec->getNameOfElement ((*plys)[k]) << std::endl;
			}
	}

	os << "#STOP" << std::endl;
	os.close ();
}

void writeAllDataToGLIFileV4 (const std::string& fname, const GeoLib::GEOObjects& geo)
{
	std::vector<std::string> geo_names;
	geo.getGeometryNames (geo_names);

	// extract path for reading external files
	const std::string path = BaseLib::extractPath(fname);

	std::ofstream os (fname.c_str());

	size_t pnts_offset (0);
	std::vector<size_t> pnts_id_offset;
	pnts_id_offset.push_back (0);

	// writing all points
	os << "#POINTS" << std::endl;
	for (size_t j(0); j < geo_names.size(); j++)
	{
		os.precision (20);
		GeoLib::PointVec const* const pnt_vec(geo.getPointVecObj(geo_names[j]));
		std::vector<GeoLib::Point*> const* const pnts (pnt_vec->getVector());
		if (pnts) {
			std::string pnt_name;
			const size_t n_pnts(pnts->size());
			for (size_t k(0); k < n_pnts; k++) {
				os << pnts_offset + k << " " << *((*pnts)[k]) << std::flush;
				if (pnt_vec->getNameOfElementByID(k, pnt_name)) {
					os << " $NAME " << pnt_name << std::flush;
				}
				os << std::endl;
			}
			pnts_offset += pnts->size();
			pnts_id_offset.push_back (pnts_offset);
		}
	}

	INFO("GeoLib::writeAllDataToGLIFileV4(): wrote %d points.", pnts_offset);

	// writing all stations
	std::vector<std::string> stn_names;
	geo.getStationVectorNames (stn_names);
	for (size_t j(0); j < stn_names.size(); j++)
	{
		os.precision (20);
		const std::vector<GeoLib::Point*>* pnts (geo.getStationVec(stn_names[j]));
		if (pnts)
		{
			for (size_t k(0); k < pnts->size(); k++)
				os << k + pnts_offset << " " << *((*pnts)[k]) << " $NAME " <<
				static_cast<GeoLib::Station*>((*pnts)[k])->getName() <<
				std::endl;
			pnts_offset += pnts->size();
			pnts_id_offset.push_back (pnts_offset);
		}
	}

	size_t plys_cnt (0);

	// writing all polylines
	for (size_t j(0); j < geo_names.size(); j++)
	{
		const GeoLib::PolylineVec* plys_vec (geo.getPolylineVecObj (geo_names[j]));
		if (plys_vec) {
			const std::vector<GeoLib::Polyline*>* plys (plys_vec->getVector());
			for (size_t k(0); k < plys->size(); k++) {
				os << "#POLYLINE" << std::endl;
				std::string ply_name;
				if (plys_vec->getNameOfElementByID (plys_cnt, ply_name))
					os << "\t$NAME " << std::endl << "\t\t" << ply_name <<
					std::endl;
				else
					os << "\t$NAME " << std::endl << "\t\t" << geo_names[j] <<
					"-" << plys_cnt << std::endl;
				os << "\t$POINTS" << std::endl;
				for (size_t l(0); l < (*plys)[k]->getNumberOfPoints(); l++)
					os << "\t\t" << pnts_id_offset[j] +
					((*plys)[k])->getPointID(l) << std::endl;
				plys_cnt++;
			}
		}
	}

	// writing surfaces as TIN files
	size_t sfcs_cnt (0);
	for (size_t j(0); j < geo_names.size(); j++)
	{
		const GeoLib::SurfaceVec* sfcs_vec (geo.getSurfaceVecObj (geo_names[j]));
			if (sfcs_vec) {
				const std::vector<GeoLib::Surface*>* sfcs (sfcs_vec->getVector());
				for (size_t k(0); k < sfcs->size(); k++)
				{
					os << "#SURFACE" << std::endl;
					std::string sfc_name(path);
					if (sfcs_vec->getNameOfElementByID (sfcs_cnt, sfc_name)) {
						os << "\t$NAME " << std::endl << "\t\t" << sfc_name << std::endl;
					} else {
						os << "\t$NAME " << std::endl << "\t\t" << sfcs_cnt << std::endl;
					sfc_name += BaseLib::number2str (sfcs_cnt);
				}
				sfc_name += ".tin";
				os << "\t$TIN" << std::endl;
				os << "\t\t" << sfc_name << std::endl;
				// create tin file
				std::ofstream tin_os (sfc_name.c_str());
				GeoLib::Surface const& sfc (*(*sfcs)[k]);
				const size_t n_tris (sfc.getNTriangles());
					for (size_t l(0); l < n_tris; l++) {
						GeoLib::Triangle const& tri (*(sfc[l]));
					tin_os << l << " " << *(tri.getPoint(0)) << " " <<
					*(tri.getPoint(1)) << " " << *(tri.getPoint(2)) <<
					std::endl;
				}
				tin_os.close();

				sfcs_cnt++;
			}
		}
	}

	os << "#STOP" << std::endl;
	os.close ();
}
} // end namespace
