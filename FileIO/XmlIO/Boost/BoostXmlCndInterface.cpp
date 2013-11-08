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
#include <fstream>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>

#include "logog/include/logog.hpp"

#include "BoostXmlCndInterface.h"

// BaseLib
#include "StringTools.h"

// OGS
#include "BoundaryCondition.h"

namespace FileIO
{

BoostXmlCndInterface::BoostXmlCndInterface(ProjectData & project_data) :
		_project_data(project_data)
{}

bool BoostXmlCndInterface::readFile(const std::string &fname)
{
	std::ifstream in(fname.c_str());
	if (!in) {
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
			boundary_condition_nodes) {
		if (boundary_condition_node.first.compare("BC") != 0)
			continue;

		// parse attribute of boundary condition
		std::string const& geometry_name = boundary_condition_node.second.get<std::string>("<xmlattr>.geometry");

		if (_project_data.getGEOObjects()->exists(geometry_name) == -1) {
			ERR("BoostXmlCndInterface::readBoundaryConditions(): Associated geometry \"%s\" not found.",
				geometry_name.c_str());
			return;
		}
		// create instance
		BoundaryCondition *bc(new BoundaryCondition(geometry_name));

		// parse tags of boundary condition
		BOOST_FOREACH(ptree::value_type const & boundary_condition_tag, boundary_condition_node.second) {
			if (boundary_condition_tag.first.compare("Process") == 0) {
				std::string pcs_type, primary_variable;
				readProcessTag(boundary_condition_tag.second, pcs_type, primary_variable);
				bc->setProcessType(FiniteElement::convertProcessType(pcs_type));
				bc->setProcessPrimaryVariable(FiniteElement::convertPrimaryVariable(primary_variable));
			}
			if (boundary_condition_tag.first.compare("Geometry") == 0) {
				std::string geo_type, geo_name;
				readGeometryTag(boundary_condition_tag.second, geo_type, geo_name);
				bc->setGeoName(geo_name);
				bc->setGeoType(GeoLib::convertGeoType(geo_type));
			}
			if (boundary_condition_tag.first.compare("Distribution") == 0) {
				readDistributionTag(boundary_condition_tag.second, bc);
			}
		}
		_project_data.addCondition(bc);
	}
}

void BoostXmlCndInterface::readProcessTag(boost::property_tree::ptree const& pcs_tags,
		std::string &pcs_type, std::string &primary_variable) const
{
	using boost::property_tree::ptree;
	BOOST_FOREACH(ptree::value_type const & pcs_tag, pcs_tags) {
		if (pcs_tag.first.compare("Type") == 0) {
			pcs_type = pcs_tag.second.data();
		}
		if (pcs_tag.first.compare("Variable") == 0) {
			primary_variable = pcs_tag.second.data();
		}
	}
}

void BoostXmlCndInterface::readGeometryTag(boost::property_tree::ptree const& geometry_tags,
		std::string &geo_type, std::string &geo_name) const
{
	using boost::property_tree::ptree;
	BOOST_FOREACH(ptree::value_type const & geo_tag, geometry_tags) {
		if (geo_tag.first.compare("Type") == 0) {
			geo_type = geo_tag.second.data();
		}
		if (geo_tag.first.compare("Name") == 0) {
			geo_name = geo_tag.second.data();
		}
	}
}

void BoostXmlCndInterface::readDistributionTag(boost::property_tree::ptree const& distribution_tags,
		FEMCondition * cond) const
{
	using boost::property_tree::ptree;
	BOOST_FOREACH(ptree::value_type const & dis_tag, distribution_tags) {

		if (dis_tag.first.compare("Type") == 0) {
			cond->setProcessDistributionType(
					FiniteElement::convertDisType(dis_tag.second.data()));
		}

		if (dis_tag.first.compare("Value") == 0) {
			FiniteElement::DistributionType const& dt(cond->getProcessDistributionType());

			if (dt == FiniteElement::CONSTANT || dt == FiniteElement::CONSTANT_NEUMANN) {
				cond->setConstantDisValue(BaseLib::str2number<double>(dis_tag.second.data()));
				return;
			}

			if (dt == FiniteElement::LINEAR || dt == FiniteElement::LINEAR_NEUMANN
					|| dt == FiniteElement::DIRECT) {
				std::vector<std::size_t> dis_node_ids;
				std::vector<double> dis_values;

				boost::tokenizer<> tok(dis_tag.second.data());
				for (boost::tokenizer<>::iterator tok_it=tok.begin(); tok_it!=tok.end(); ) {
					dis_node_ids.push_back(BaseLib::str2number<std::size_t>(*tok_it));
					tok_it++;
					dis_values.push_back(BaseLib::str2number<double>(*tok_it));
					tok_it++;
				}
				cond->setDisValues(dis_node_ids, dis_values);
				return;
			}

			ERR("BoostXmlCndInterface::readDistributionTag(): Distribution type not supported.")
		}
	}
}

int BoostXmlCndInterface::write(std::ostream& stream)
{
	stream << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"; // xml definition
	stream << "<?xml-stylesheet type=\"text/xsl\" href=\"OpenGeoSysCND.xsl\"?>\n\n"; // stylefile definition

	// create a DOM tree for writing it to file
	using boost::property_tree::ptree;
	ptree pt;

	pt.add("OpenGeoSysCond", "OpenGeoSysCond");
	pt.put("<xmlattr>.xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
	pt.put("<xmlattr>.xsi:noNamespaceSchemaLocation", "http://www.opengeosys.org/images/xsd/OpenGeoSysCND.xsd");
	pt.put("<xmlattr>.xmlns:ogs", "http://www.opengeosys.net");

	write_xml(stream, pt);

	return 0;
}

} // end namespace FileIO
