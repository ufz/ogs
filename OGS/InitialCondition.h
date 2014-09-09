/**
 * \copyright
 * Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef OGS_INITIAL_CONDITION_H_
#define OGS_INITIAL_CONDITION_H_

#include <boost/property_tree/ptree.hpp>
#include "logog/include/logog.hpp"

#include "GeoLib/GEOObjects.h"
#include "MeshLib/Mesh.h"

namespace OGS
{

class InitialCondition
{
public:
    virtual ~InitialCondition() = default;
};


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

private:
    double _value;
};


}   // namespace OGS

#endif  // OGS_INITIAL_CONDITION_H_
