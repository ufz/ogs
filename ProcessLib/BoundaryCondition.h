/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_BOUNDARY_CONDITION_H_
#define PROCESS_LIB_BOUNDARY_CONDITION_H_

#include <boost/property_tree/ptree.hpp>
#include "logog/include/logog.hpp"

namespace GeoLib
{
    class GeoObject;
}

namespace ProcessLib
{

class BoundaryCondition
{
public:
    BoundaryCondition(GeoLib::GeoObject const* const geometry)
        : _geometry(geometry)
    { }

    virtual ~BoundaryCondition() = default;

private:
    GeoLib::GeoObject const* const _geometry;
};

/// The UniformDirichletBoundaryCondition class describes a constant in space
/// and time Dirichlet boundary condition.
/// The expected parameter in the passed configuration is "value" which, when
/// not present defaults to zero.
class UniformDirichletBoundaryCondition : public BoundaryCondition
{
    using ConfigTree = boost::property_tree::ptree;
public:
    UniformDirichletBoundaryCondition(GeoLib::GeoObject const* const geometry,
            ConfigTree const& config)
        : BoundaryCondition(geometry)
    {
        DBUG("Constructing UniformDirichletBoundaryCondition from config.");

        _value = config.get<double>("value", 0);
        DBUG("Using value %g", _value);
    }

private:
    double _value;
};


}   // namespace ProcessLib

#endif  // PROCESS_LIB_BOUNDARY_CONDITION_H_
