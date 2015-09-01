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

}  // namespace ProcessLib
