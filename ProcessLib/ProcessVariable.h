/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_PROCESS_VARIABLE_H_
#define PROCESS_LIB_PROCESS_VARIABLE_H_

#include <boost/property_tree/ptree.hpp>

#include "GeoLib/GEOObjects.h"
#include "MeshLib/Mesh.h"

namespace ProcessLib
{
    class BoundaryCondition;
    class InitialCondition;
}

namespace ProcessLib
{

/// A named process variable. Its properties includes the mesh, and the initial
/// and boundary conditions.
class ProcessVariable
{
    using ConfigTree = boost::property_tree::ptree;
public:
    ProcessVariable(ConfigTree const& config, MeshLib::Mesh const& mesh,
            GeoLib::GEOObjects const& geometries);

    ~ProcessVariable();

    std::string const& getName() const;

    /// Returns a mesh on which the process variable is defined.
    MeshLib::Mesh const& getMesh() const;

    /// Const iterator over boundary conditions of the process variable.
    using BoundaryConditionCI = std::vector<BoundaryCondition*>::const_iterator;

    /// Returns a BoundaryConditionCI iterator to the beginning.
    BoundaryConditionCI
    beginBoundaryConditions() const
    {
        return _boundary_conditions.cbegin();
    }

    /// Returns a past-the-end BoundaryConditionCI iterator.
    BoundaryConditionCI
    endBoundaryConditions() const
    {
        return _boundary_conditions.cend();
    }

private:
    std::string const _name;
    MeshLib::Mesh const& _mesh;
    InitialCondition* _initial_condition;
    std::vector<BoundaryCondition*> _boundary_conditions;
};

}   // namespace ProcessLib

#endif  // PROCESS_LIB_PROCESS_VARIABLE_H_
