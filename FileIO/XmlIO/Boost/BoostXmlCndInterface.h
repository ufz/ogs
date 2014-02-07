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

namespace FileIO
{

/**
 * \brief Reads and writes FEM Conditions to and from XML files using the boost XML parser.
 */
class BoostXmlCndInterface : public XMLInterface
{
public:
	explicit BoostXmlCndInterface(ProjectData & project);
	virtual ~BoostXmlCndInterface()	{}

	/// Reads an xml-file containing FEM Conditions such as Boundary- or Initial Conditions
	bool readFile(const std::string &fname);

	void setConditionType(FEMCondition::CondType type) { _type = type; }

protected:
	/// @return true on success, else false
	bool write(std::ostream& stream);

private:
	/// Read the details of a boundary condition from an xml-file
	void readBoundaryConditions(boost::property_tree::ptree const& boundary_condition_nodes);

	/// Read details on process parameters
	void readProcessTag(boost::property_tree::ptree const& pcs_tags,
			std::string &pcs_type, std::string &primary_variable) const;

	/// Read details on geometric parameters
	void readGeometryTag(boost::property_tree::ptree const& geometry_tags,
			std::string &geo_type, std::string &geo_name) const;

	/// Read details on distribution parameters
	void readDistributionTag(boost::property_tree::ptree const& distribution_tags,
			FEMCondition * bc) const;

	FEMCondition::CondType _type;
	ProjectData & _project_data;
};

} // end namespace FileIO

#endif /* BOOSTXMLCNDINTERFACE_H_ */
