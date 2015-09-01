/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_INITIAL_CONDITION_H_
#define PROCESS_LIB_INITIAL_CONDITION_H_

#include <boost/property_tree/ptree.hpp>
#include "logog/include/logog.hpp"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

namespace ProcessLib
{
/// The InitialCondition is a base class for spatial distributions of values
/// defined on mesh nodes.
class InitialCondition
{
public:
	virtual ~InitialCondition() = default;
	virtual double getValue(MeshLib::Node const&) const = 0;
};

/// Uniform value initial condition
class UniformInitialCondition : public InitialCondition
{
	using ConfigTree = boost::property_tree::ptree;

public:
	UniformInitialCondition(ConfigTree const& config)
	{
		DBUG("Constructing Uniform initial condition");

		_value = config.get<double>("value", 0);
		DBUG("Read value %g", _value);
	}

	virtual double getValue(MeshLib::Node const&) const override
	{
		return _value;
	}

private:
	double _value;
};

}  // namespace ProcessLib

#endif  // PROCESS_LIB_INITIAL_CONDITION_H_
