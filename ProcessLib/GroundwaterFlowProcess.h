/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_
#define PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_

#include <memory>

#include <boost/property_tree/ptree.hpp>
#include "logog/include/logog.hpp"

#include "AssemblerLib/LocalDataInitializer.h"
#include "AssemblerLib/LocalToGlobalIndexMap.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSubset.h"
#include "MeshLib/MeshSubsets.h"

#include "GroundwaterFlowFEM.h"
#include "ProcessVariable.h"

namespace ProcessLib
{

template<typename GlobalSetup>
class GroundwaterFlowProcess : public Process
{
    using ConfigTree = boost::property_tree::ptree;

public:
    GroundwaterFlowProcess(MeshLib::Mesh const& mesh,
            std::vector<ProcessVariable> const& variables,
            ConfigTree const& config)
        : Process(mesh, 1)
    {
        DBUG("Create GroundwaterFlowProcess.");

        // Find the corresponding process variable.
        std::string const name = config.get<std::string>("process_variable");

        auto const& variable = std::find_if(variables.cbegin(), variables.cend(),
                [&name](ProcessVariable const& v) {
                    return v.getName() == name;
                });

        if (variable == variables.end())
            ERR("Expected process variable \'%s\' not found in provided variables list.",
                name.c_str());

        DBUG("Associate hydraulic_head with process variable \'%s\'.",
            name.c_str());
        _hydraulic_head = &*variable;

    }

    void initialize()
    {
        DBUG("Construct dof mappings.");
        // Create single component dof in every of the mesh's nodes.
        _mesh_subset_all_nodes = new MeshLib::MeshSubset(_mesh, _mesh.getNodes());

        // Collect the mesh subsets in a vector.
        _all_mesh_subsets.push_back(new MeshLib::MeshSubsets(_mesh_subset_all_nodes));

        _local_to_global_index_map.reset(
            new AssemblerLib::LocalToGlobalIndexMap(_all_mesh_subsets));

        //DBUG("Create global assembler.");
        //_global_assembler.reset(
        //    new GlobalAssembler(*_A, *_rhs, *_local_to_global_index_map));
    }

    void solve()
    {
    }

    void post()
    {
    }

    ~GroundwaterFlowProcess()
    {
        for (auto p : _all_mesh_subsets)
            delete p;

        delete _mesh_subset_all_nodes;
    }

private:
    ProcessVariable const* _hydraulic_head = nullptr;

    MeshLib::MeshSubset const* _mesh_subset_all_nodes = nullptr;
    std::vector<MeshLib::MeshSubsets*> _all_mesh_subsets;

    GlobalSetup _global_setup;

    std::unique_ptr<AssemblerLib::LocalToGlobalIndexMap> _local_to_global_index_map;
};

}   // namespace ProcessLib

#endif  // PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_
