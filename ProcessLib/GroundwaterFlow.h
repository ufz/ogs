/**
 * \copyright
 * Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_GROUNDWATERFLOW_H_
#define PROCESS_LIB_GROUNDWATERFLOW_H_

#include <boost/property_tree/ptree.hpp>
#include "logog/include/logog.hpp"

#include "AssemblerLib/MeshComponentMap.h"
#include "GeoLib/GEOObjects.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSubset.h"
#include "MeshLib/MeshSubsets.h"
#include "OGS/ProcessVariable.h"

namespace ProcessLib
{

class GroundwaterFlowProcess : public Process
{
    using ConfigTree = boost::property_tree::ptree;

public:
    GroundwaterFlowProcess(MeshLib::Mesh const& mesh,
            std::vector<OGS::ProcessVariable> const& variables,
            ConfigTree const& config)
        : Process(mesh, 1)
    {
        DBUG("Create GroundwaterFlowProcess.");

        // Find the corresponding process variable.
        std::string const name = config.get<std::string>("process_variable");

        auto const& variable = std::find_if(variables.cbegin(), variables.cend(),
                [&name](OGS::ProcessVariable const& v) {
                    return v.getName() == name;
                });

        if (variable == variables.end())
            ERR("Expected field \'%s\' not found in provided fields.",
                name.c_str());

        DBUG("Associate hydraulic_head with field named \'%s\'.",
            name.c_str());
        _hydraulic_head = &*variable;

    }

    void initialize()
    {
        // Create mesh's subset using all nodes of the mesh.
        _mesh_subset_all_nodes = new MeshLib::MeshSubset(_mesh, _mesh.getNodes());

        // Define a mesh item composition in a vector.
        _all_mesh_subsets.push_back(new MeshLib::MeshSubsets(_mesh_subset_all_nodes));
        AssemblerLib::MeshComponentMap mesh_component_map(_all_mesh_subsets,
                AssemblerLib::ComponentOrder::BY_COMPONENT);
    }

    ~GroundwaterFlowProcess()
    {
        for (auto p : _all_mesh_subsets)
            delete p;

        delete _mesh_subset_all_nodes;
    }


private:
    OGS::ProcessVariable const* _hydraulic_head = nullptr;

    MeshLib::MeshSubset const* _mesh_subset_all_nodes = nullptr;
    std::vector<MeshLib::MeshSubsets*> _all_mesh_subsets;
};

}   // namespace ProcessLib

#endif  // PROCESS_LIB_GROUNDWATERFLOW_H_
