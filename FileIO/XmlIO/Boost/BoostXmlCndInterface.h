/**
 * @file BoostXmlCndInterface.h
 * @author git blame BoostXmlCndInterface.h
 * @date Oct 14, 2013
 * @brief Class BoostXmlCndInterface is for reading FEM conditions
 * 	(initial, boundary conditions or source terms). Implementation uses boost.
 *
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */
#ifndef BOOSTXMLCNDINTERFACE_H_
#define BOOSTXMLCNDINTERFACE_H_

#include <string>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/optional.hpp>

#include "../XMLInterface.h"

class FEMCondition;
class Projectdata;

typedef boost::optional<boost::property_tree::ptree const&> OptionalPtree;

namespace FileIO
{

class BoostXmlCndInterface : public XMLInterface
{
public:
	BoostXmlCndInterface(ProjectData* project);
	virtual ~BoostXmlCndInterface()	{}

	/// Reads an xml-file containing FEM Conditions such as Boundary- or Initial Conditions
	bool readFile(const std::string &fname);

protected:
	int write(std::ostream& stream);

private:
	void readBoundaryConditions(boost::property_tree::ptree const& boundary_condition_nodes);
	void readProcessTag(boost::property_tree::ptree const& pcs_tags,
			std::string &pcs_type, std::string &primary_variable) const;
	void readGeometryTag(boost::property_tree::ptree const& geometry_tags,
			std::string &geo_type, std::string &geo_name) const;
	void readDistributionTag(boost::property_tree::ptree const& distribution_tags,
			FEMCondition * bc) const;

	ProjectData* _project_data;
};

} // end namespace FileIO

#endif /* BOOSTXMLCNDINTERFACE_H_ */
