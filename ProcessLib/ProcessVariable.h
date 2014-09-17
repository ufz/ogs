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

private:
    std::string const _name;
    MeshLib::Mesh const& _mesh;
    InitialCondition* _initial_condition;
    std::vector<BoundaryCondition*> _boundary_conditions;
};

}   // namespace ProcessLib

#endif  // PROCESS_LIB_PROCESS_VARIABLE_H_
