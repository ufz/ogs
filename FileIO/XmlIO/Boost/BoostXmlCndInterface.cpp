/**
 * @file BoostXmlCndInterface.cpp
 * @author git blame BoostXmlCndInterface.cpp
 * @date Oct 14, 2013
 * @brief 
 *
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */
#include "BoostXmlCndInterface.h"

#include <fstream>

#include <boost/foreach.hpp>

#include "logog/include/logog.hpp"

#include "BoundaryCondition.h"

namespace FileIO
{

BoostXmlCndInterface::BoostXmlCndInterface(ProjectData* project_data) :
		_project_data(project_data)
{}

bool BoostXmlCndInterface::readFile(const std::string &fname)
{
	std::ifstream in(fname.c_str());
	if (in.fail()) {
		ERR("BoostXmlCndInterface::readFile(): Can't open xml-file %s.", fname.c_str());
		return false;
	}

	// build DOM tree
	using boost::property_tree::ptree;
	ptree pt;
	read_xml(in, pt);

	ptree const& root_node = pt.get_child("OpenGeoSysCond");

	BOOST_FOREACH(ptree::value_type const & conditions_type, root_node) {
		if (conditions_type.first.compare("BoundaryConditions") == 0) {
			readBoundaryConditions(conditions_type.second);
		}
	}

	return true;
}

void BoostXmlCndInterface::readBoundaryConditions(
		boost::property_tree::ptree const& boundary_condition_nodes)
{
	using boost::property_tree::ptree;
	BOOST_FOREACH(ptree::value_type const & boundary_condition_node,
			boundary_condition_nodes)
	{
//		if (boundary_condition_node.first.compare("BC") == 0) {
//			// parse attribute of boundary condition
//			std::string const& geometry_name = boundary_condition_node.get<std::string>("<xmlattr>.geometry");
//
//			// create instance
//			BoundaryCondition *bc(new BoundaryCondition(geometry_name));
//			// ToDo: check if geometry exists
//
//			// parse tags of boundary condition
//			BOOST_FOREACH(ptree::value_type const & boundary_condition_tag, *boundary_condition_node) {
//				if (boundary_condition_tag.first.compare("Process") == 0) {
//					std::string type, primary_variable;
//					readProcessInfo(boundary_condition_tag, type, primary_variable);
//					bc->setProcessType(FiniteElement::convertProcessType(type));
//					bc->setProcessPrimaryVariable(FiniteElement::convertPrimaryVariable(primary_variable));
//				}
//				if (boundary_condition_tag.first.compare("Geometry") == 0) {
//					std::string geo_type, geo_name;
//					readGeometryInfo(boundary_condition_tag, geo_type, geo_name);
//					bc->setGeoName(geo_name);
//					bc->setGeoType(GeoLib::convertGeoType(geo_type));
//				}
//			}
//		}
	}
}

int BoostXmlCndInterface::write(std::ostream& stream)
{
	return 0;
}

} // end namespace FileIO
