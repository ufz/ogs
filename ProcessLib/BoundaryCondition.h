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

#include <algorithm>
#include <vector>

#include <boost/property_tree/ptree.hpp>
#include "logog/include/logog.hpp"

#include "MeshGeoToolsLib/MeshNodeSearcher.h"

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

    /// Initialization interface of Dirichlet type boundary conditions.
    /// Fills in global_ids of the particular geometry of the boundary condition
    /// and the corresponding values.
    /// The ids are appended to the global_ids and also the values.
    virtual void initialize(
            MeshGeoToolsLib::MeshNodeSearcher& searcher,
            std::vector<std::size_t>& global_ids, std::vector<double>& values) = 0;

protected:
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

    /// Initialize Dirichlet type boundary conditions.
    /// Fills in global_ids of the particular geometry of the boundary condition
    /// and the corresponding values.
    /// The ids are appended to the global_ids and the values are filled with
    /// the constant _value.
    void initialize(MeshGeoToolsLib::MeshNodeSearcher& searcher,
            std::vector<std::size_t>& global_ids, std::vector<double>& values)
    {
        // Find nodes' ids on the given mesh on which this boundary condition
        // is defined.
        std::vector<std::size_t> ids = searcher.getMeshNodeIDs(*_geometry);

        // Append node ids.
        global_ids.reserve(global_ids.size() + ids.size());
        std::copy(ids.cbegin(), ids.cend(), std::back_inserter(global_ids));

        // Fill values.
        values.reserve(values.size() + ids.size());
        std::fill_n(std::back_inserter(values), ids.size(), _value);
    }

private:
    double _value;
};


}   // namespace ProcessLib

#endif  // PROCESS_LIB_BOUNDARY_CONDITION_H_
