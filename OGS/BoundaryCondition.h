/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef OGS_BOUNDARY_CONDITION_H_
#define OGS_BOUNDARY_CONDITION_H_

#include <boost/property_tree/ptree.hpp>
#include "logog/include/logog.hpp"

namespace GeoLib
{
    class GeoObject;
}

namespace OGS
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


class DirichletBoundaryCondition : public BoundaryCondition
{
    using ConfigTree = boost::property_tree::ptree;
public:
    DirichletBoundaryCondition(GeoLib::GeoObject const* const geometry,
            ConfigTree const& config)
        : BoundaryCondition(geometry)
    {
        DBUG("Constructing Dirichlet boundary condition");

        _value = config.get<double>("value", 0);
        DBUG("Read value %g", _value);
    }

private:
    double _value;
};


}   // namespace OGS

#endif  // OGS_BOUNDARY_CONDITION_H_
