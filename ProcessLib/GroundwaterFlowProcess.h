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

#include "AssemblerLib/LocalAssemblerBuilder.h"
#include "AssemblerLib/VectorMatrixAssembler.h"
#include "AssemblerLib/LocalDataInitializer.h"
#include "AssemblerLib/LocalToGlobalIndexMap.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSubset.h"
#include "MeshLib/MeshSubsets.h"

#include "NumLib/Fem/Integration/IntegrationGaussRegular.h"

#include "GroundwaterFlowFEM.h"
#include "ProcessVariable.h"

namespace ProcessLib
{

template<typename GlobalSetup>
class GroundwaterFlowProcess : public Process
{
    using ConfigTree = boost::property_tree::ptree;

    template <typename ShapeFunction_>
    using IntegrationPolicy = NumLib::IntegrationGaussRegular<ShapeFunction_::DIM>;
    unsigned const _integration_order = 2;

public:
    GroundwaterFlowProcess(MeshLib::Mesh const& mesh,
            std::vector<ProcessVariable> const& variables,
            ConfigTree const& config)
        : Process(mesh)
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
        DBUG("Initialize GroundwaterFlowProcess.");

        DBUG("Construct dof mappings.");
        // Create single component dof in every of the mesh's nodes.
        _mesh_subset_all_nodes = new MeshLib::MeshSubset(_mesh, _mesh.getNodes());

        // Collect the mesh subsets in a vector.
        _all_mesh_subsets.push_back(new MeshLib::MeshSubsets(_mesh_subset_all_nodes));

        _local_to_global_index_map.reset(
            new AssemblerLib::LocalToGlobalIndexMap(_all_mesh_subsets));

        DBUG("Allocate global matrix, vectors, and linear solver.");
        _A.reset(_global_setup.createMatrix(_local_to_global_index_map->dofSize()));
        _x.reset(_global_setup.createVector(_local_to_global_index_map->dofSize()));
        _rhs.reset(_global_setup.createVector(_local_to_global_index_map->dofSize()));
        _linearSolver.reset(new typename GlobalSetup::LinearSolver(*_A));

        DBUG("Create local assemblers.");
        // Populate the vector of local assemblers.
        _local_assemblers.resize(_mesh.getNElements());
        // Shape matrices initializer
        using LocalDataInitializer = AssemblerLib::LocalDataInitializer<
            GroundwaterFlow::LocalAssemblerDataInterface,
            GroundwaterFlow::LocalAssemblerData,
            IntegrationPolicy,
            typename GlobalSetup::MatrixType,
            typename GlobalSetup::VectorType>;

        LocalDataInitializer initializer;

        using LocalAssemblerBuilder =
            AssemblerLib::LocalAssemblerBuilder<
                MeshLib::Element,
                LocalDataInitializer>;

        LocalAssemblerBuilder local_asm_builder(
            initializer, *_local_to_global_index_map);

        DBUG("Calling local assembler builder for all mesh elements.");
        _global_setup.execute(
                local_asm_builder,
                _mesh.getElements(),
                _local_assemblers,
                1e-6,   // hydraulic conductivity
                _integration_order);

        DBUG("Create global assembler.");
        _global_assembler.reset(
            new GlobalAssembler(*_A, *_rhs, *_local_to_global_index_map));
    }

    void solve()
    {
        DBUG("Solve GroundwaterFlowProcess.");

        _A->setZero();
        *_rhs = 0;   // This resets the whole vector.

        // Call global assembler for each local assembly item.
        _global_setup.execute(*_global_assembler, _local_assemblers);

        _linearSolver->solve(*_rhs, *_x);
    }

    void post()
    {
        DBUG("Postprocessing GroundwaterFlowProcess.");
        // Postprocessing of the linear system of equations solver results:
        // For example, write _x to _hydraulic_head or convert to velocity etc.
    }

    ~GroundwaterFlowProcess()
    {
        for (auto p : _local_assemblers)
            delete p;

        for (auto p : _all_mesh_subsets)
            delete p;

        delete _mesh_subset_all_nodes;
    }

private:
    ProcessVariable const* _hydraulic_head = nullptr;

    MeshLib::MeshSubset const* _mesh_subset_all_nodes = nullptr;
    std::vector<MeshLib::MeshSubsets*> _all_mesh_subsets;

    GlobalSetup _global_setup;
    std::unique_ptr<typename GlobalSetup::LinearSolver> _linearSolver;
    std::unique_ptr<typename GlobalSetup::MatrixType> _A;
    std::unique_ptr<typename GlobalSetup::VectorType> _rhs;
    std::unique_ptr<typename GlobalSetup::VectorType> _x;

    std::vector<GroundwaterFlow::LocalAssemblerDataInterface<
        typename GlobalSetup::MatrixType, typename GlobalSetup::VectorType>*>
            _local_assemblers;

    using GlobalAssembler =
        AssemblerLib::VectorMatrixAssembler<
            typename GlobalSetup::MatrixType,
            typename GlobalSetup::VectorType>;


    std::unique_ptr<AssemblerLib::LocalToGlobalIndexMap> _local_to_global_index_map;

    std::unique_ptr<GlobalAssembler> _global_assembler;
};

}   // namespace ProcessLib

#endif  // PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_
