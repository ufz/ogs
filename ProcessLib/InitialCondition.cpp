/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "InitialCondition.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/optional.hpp>
#include <logog/include/logog.hpp>

#include "MathLib/Point3d.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"

namespace ProcessLib
{
using ConfigTree = boost::property_tree::ptree;

std::unique_ptr<InitialCondition> createUniformInitialCondition(
    ConfigTree const& config)
{
	auto value = config.get_optional<double>("value");
	if (!value)
	{
		ERR("Could not find required parameter value.");
		std::abort();
	}
	DBUG("Using value %g", *value);

	return std::unique_ptr<InitialCondition>(
	    new UniformInitialCondition(*value));
}

std::unique_ptr<InitialCondition> createMeshPropertyInitialCondition(
    ConfigTree const& config, MeshLib::Mesh const& mesh)
{
	auto field_name = config.get_optional<std::string>("field_name");
	if (!field_name)
	{
		ERR("Could not find required parameter field_name.");
		std::abort();
	}
	DBUG("Using field_name %s", field_name->c_str());

	if (!mesh.getProperties().hasPropertyVector(*field_name))
	{
		ERR("The required property %s does not exists in the mesh.",
		    field_name->c_str());
		std::abort();
	}
	auto const& property =
	    mesh.getProperties().template getPropertyVector<double>(*field_name);
	if (!property)
	{
		ERR("The required property %s is not of the requested type.",
		    field_name->c_str());
		std::abort();
	}

	return std::unique_ptr<InitialCondition>(
	    new MeshPropertyInitialCondition(*property));
}

}  // namespace ProcessLib
