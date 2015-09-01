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

#include <boost/property_tree/ptree_fwd.hpp>
#include "MeshLib/Node.h"
#include "MeshLib/PropertyVector.h"

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
public:
	UniformInitialCondition(double const value) : _value(value)
	{
	}

	virtual double getValue(MeshLib::Node const&) const override
	{
		return _value;
	}

private:
	double _value;
};

using ConfigTree = boost::property_tree::ptree;
/// Construct a UniformInitialCondition from configuration.
std::unique_ptr<InitialCondition> createUniformInitialCondition(
    ConfigTree const& config);

}  // namespace ProcessLib

#endif  // PROCESS_LIB_INITIAL_CONDITION_H_
