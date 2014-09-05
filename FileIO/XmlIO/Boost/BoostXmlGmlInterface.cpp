/**
 * \file   BoostXmlGmlInterface.cpp
 * \author Karsten Rink
 * \date   2014-01-31
 * \brief  Implementation of the BoostXmlGmlInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>

#include "BoostXmlGmlInterface.h"

#include <fstream>

#include "logog/include/logog.hpp"

#include "ProjectData.h"
#include "GEOObjects.h"


namespace FileIO
{

BoostXmlGmlInterface::BoostXmlGmlInterface(GeoLib::GEOObjects& geo_objs) :
		_geo_objects(geo_objs)
{}

bool BoostXmlGmlInterface::readFile(const std::string &fname)
{
	std::ifstream in(fname.c_str());
	if (!in) {
		ERR("BoostXmlGmlInterface::readFile(): Can't open xml-file %s.", fname.c_str());
		return false;
	}

	std::string geo_name("[NN]");

	std::vector<GeoLib::Point*>* points = new std::vector<GeoLib::Point*>;
	std::vector<GeoLib::Polyline*>* polylines = new std::vector<GeoLib::Polyline*>;
	std::vector<GeoLib::Surface*>* surfaces = new std::vector<GeoLib::Surface*>;

	std::map<std::string, std::size_t>* pnt_names = new std::map<std::string, std::size_t>;
	std::map<std::string, std::size_t>* ply_names = new std::map<std::string, std::size_t>;
	std::map<std::string, std::size_t>* sfc_names  = new std::map<std::string, std::size_t>;

	GeoLib::GEOObjects* geo_objects (&_geo_objects);


	// build DOM tree
	using boost::property_tree::ptree;
	ptree doc;
	read_xml(in, doc);

	
	if (!isGmlFile(doc))
		return false;

	ptree const & root_node = doc.get_child("OpenGeoSysGLI");
	BOOST_FOREACH( ptree::value_type const & node, root_node )
	{
		if (node.first.compare("name") == 0)
		{
			if (node.second.data().empty())
			{
				ERR("BoostXmlGmlInterface::readFile(): <name>-tag is empty.")
				return false;
			}
			geo_name = node.second.data();
		}
		else if (node.first.compare("points") == 0)
		{
			readPoints(node.second, points, pnt_names);
			geo_objects->addPointVec(points, geo_name, pnt_names);
		}
		else if (node.first.compare("polylines") == 0)
			readPolylines(node.second, polylines, points, geo_objects->getPointVecObj(geo_name)->getIDMap(), ply_names);
		else if (node.first.compare("surfaces") == 0)
			readSurfaces(node.second, surfaces, points, geo_objects->getPointVecObj(geo_name)->getIDMap(), sfc_names);
	}

	if (!polylines->empty())
	{
		geo_objects->addPolylineVec(polylines, geo_name, ply_names);
	}
	else
	{
		delete polylines;
		delete ply_names;
	}

	if (!surfaces->empty())
	{
		geo_objects->addSurfaceVec(surfaces, geo_name, sfc_names);
	}
	else
	{
		delete surfaces;
		delete sfc_names;
	}

	return true;
}

void BoostXmlGmlInterface::readPoints(boost::property_tree::ptree const & pointsRoot,
	                                  std::vector<GeoLib::Point*>* points,
	                                  std::map<std::string, std::size_t>* &pnt_names )
{
	using boost::property_tree::ptree;
	BOOST_FOREACH( ptree::value_type const & point, pointsRoot )
	{
		if (point.first.compare("point") != 0)
			continue;

		unsigned p_id = static_cast<unsigned>(point.second.get("<xmlattr>.id", std::numeric_limits<unsigned>::max()));
		double   p_x  = static_cast<double>  (point.second.get("<xmlattr>.x",  std::numeric_limits<double>::max()));
		double   p_y  = static_cast<double>  (point.second.get("<xmlattr>.y",  std::numeric_limits<double>::max()));
		double   p_z  = static_cast<double>  (point.second.get("<xmlattr>.z",  std::numeric_limits<double>::max()));
		std::string p_name = point.second.get("<xmlattr>.name", "");

		if ( p_id == std::numeric_limits<unsigned>::max() || p_x == std::numeric_limits<double>::max() || 
			 p_y  == std::numeric_limits<double>::max()   || p_z == std::numeric_limits<double>::max() ) 
			WARN("BoostXmlGmlInterface::readPoints(): Attribute missing in <point> tag. Skipping point...")
		else
		{
			_idx_map.insert (std::pair<std::size_t, std::size_t>(p_id, points->size()));
			GeoLib::Point* p = new GeoLib::Point(p_x, p_y, p_z);
			if (!p_name.empty())
				pnt_names->insert( std::pair<std::string, std::size_t>(p_name, points->size()) );
			points->push_back(p);
		}
	}

	// if names-map is empty, set it to nullptr because it is not needed
	if (pnt_names->empty())
	{
		delete pnt_names;
		pnt_names = nullptr; 
	}
}


void BoostXmlGmlInterface::readPolylines(boost::property_tree::ptree const& polylinesRoot,
	                                     std::vector<GeoLib::Polyline*>* polylines,
	                                     std::vector<GeoLib::Point*>* points,
	                                     const std::vector<std::size_t> &pnt_id_map,
	                                     std::map<std::string, std::size_t>* &ply_names )
{
	using boost::property_tree::ptree;
	BOOST_FOREACH( ptree::value_type const & polyline, polylinesRoot )
	{
		if (polyline.first.compare("polyline") != 0)
			continue;

		if (static_cast<unsigned>(polyline.second.get("<xmlattr>.id", std::numeric_limits<unsigned>::max()) == std::numeric_limits<unsigned>::max()))
			WARN("BoostXmlGmlInterface::readPolylines(): Attribute \"id\" missing in <polyline> tag. Skipping polyline...")
		else
		{
			polylines->push_back(new GeoLib::Polyline(*points));

			std::string p_name = polyline.second.get("<xmlattr>.name", "");
			if (!p_name.empty())
				ply_names->insert(std::pair<std::string, std::size_t>(p_name, polylines->size()-1));

			BOOST_FOREACH( ptree::value_type const & pnt, polyline.second )
			{
				if (pnt.first.compare("pnt") == 0)
					polylines->back()->addPoint(pnt_id_map[_idx_map[atoi(pnt.second.data().c_str())]]);
			}
		}
	}

	// if names-map is empty, set it to nullptr because it is not needed
	if (ply_names->empty())
	{
		delete ply_names;
		ply_names = nullptr; 
	}
}

void BoostXmlGmlInterface::readSurfaces(boost::property_tree::ptree const& surfacesRoot,
	                                    std::vector<GeoLib::Surface*>* surfaces,
	                                    std::vector<GeoLib::Point*>* points,
	                                    const std::vector<std::size_t> &pnt_id_map,
	                                    std::map<std::string, std::size_t>* &sfc_names )
{
	using boost::property_tree::ptree;
	BOOST_FOREACH( ptree::value_type const & surface, surfacesRoot )
	{
		if (surface.first.compare("surface") != 0)
			continue;

		if (static_cast<unsigned>(surface.second.get("<xmlattr>.id", std::numeric_limits<unsigned>::max()) == std::numeric_limits<unsigned>::max()))
		{
			WARN("BoostXmlGmlInterface::readSurfaces(): Attribute \"id\" missing in <surface> tag. Skipping surface...")
			continue;
		}

		surfaces->push_back(new GeoLib::Surface(*points));

		std::string s_name = surface.second.get("<xmlattr>.name", "");
		if (!s_name.empty())
			sfc_names->insert(std::pair<std::string, std::size_t>(s_name, surfaces->size()-1));

		BOOST_FOREACH( ptree::value_type const & element, surface.second )
		{
			if (element.first.compare("element") != 0)
				continue;

			unsigned p1_attr = static_cast<unsigned>(element.second.get("<xmlattr>.p1", std::numeric_limits<unsigned>::max()));
			unsigned p2_attr = static_cast<unsigned>(element.second.get("<xmlattr>.p2", std::numeric_limits<unsigned>::max()));
			unsigned p3_attr = static_cast<unsigned>(element.second.get("<xmlattr>.p3", std::numeric_limits<unsigned>::max()));

			if (p1_attr == std::numeric_limits<unsigned>::max() || p2_attr == std::numeric_limits<unsigned>::max() || p3_attr == std::numeric_limits<unsigned>::max())
				WARN("BoostXmlGmlInterface::readSurfaces(): Attribute missing in <element> tag. Skipping triangle...");
			{
				std::size_t p1 = pnt_id_map[_idx_map[p1_attr]];
				std::size_t p2 = pnt_id_map[_idx_map[p2_attr]];
				std::size_t p3 = pnt_id_map[_idx_map[p3_attr]];
				surfaces->back()->addTriangle(p1,p2,p3);
			}
		}
	}

	// if names-map is empty, set it to nullptr because it is not needed
	if (sfc_names->empty())
	{
		delete sfc_names;
		sfc_names = nullptr; 
	}
}

bool BoostXmlGmlInterface::isGmlFile(const boost::property_tree::ptree &root) const
{
	if (!root.get_child_optional("OpenGeoSysGLI"))
	{
		ERR("BoostXmlGmlInterface::isGmlFile(): Not a GML file.");
		return false;
	}
	return true;
}

bool BoostXmlGmlInterface::write()
{
	INFO ("Writing XML geometry is not implemented here. Please use the Qt XML class for this functionality.");
	return false;
}

} // end namespace FileIO
