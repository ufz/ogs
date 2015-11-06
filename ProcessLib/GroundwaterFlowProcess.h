/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_
#define PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_

#include <cassert>
#include <memory>

#include <boost/optional.hpp>
#include <boost/algorithm/string/erase.hpp>

#include "logog/include/logog.hpp"
#include "BaseLib/ConfigTree.h"

#ifdef USE_PETSC
#include "MeshLib/NodePartitionedMesh.h"
#include "MathLib/LinAlg/PETSc/PETScMatrixOption.h"
#endif

#include "AssemblerLib/LocalAssemblerBuilder.h"
#include "AssemblerLib/VectorMatrixAssembler.h"
#include "AssemblerLib/LocalDataInitializer.h"
#include "AssemblerLib/LocalToGlobalIndexMap.h"

#include "FileIO/VtkIO/VtuInterface.h"

#include "MathLib/LinAlg/ApplyKnownSolution.h"
#include "MathLib/LinAlg/SetMatrixSparsity.h"

#include "MeshLib/MeshSubset.h"
#include "MeshLib/MeshSubsets.h"
#include "MeshLib/NodeAdjacencyTable.h"
#include "MeshGeoToolsLib/MeshNodeSearcher.h"

#include "UniformDirichletBoundaryCondition.h"

#include "GroundwaterFlowFEM.h"
#include "NeumannBcAssembler.h"
#include "NeumannBc.h"
#include "Parameter.h"
#include "Process.h"
#include "ProcessVariable.h"

namespace MeshLib
{
    class Element;
    class Mesh;
    template <typename PROP_VAL_TYPE> class PropertyVector;
};

namespace ProcessLib
{

template<typename GlobalSetup>
class GroundwaterFlowProcess : public Process
{
    unsigned const _integration_order = 2;

public:
    GroundwaterFlowProcess(MeshLib::Mesh& mesh,
            std::vector<ProcessVariable> const& variables,
            std::vector<std::unique_ptr<ParameterBase>> const& parameters,
            BaseLib::ConfigTree const& config)
        : Process(mesh)
    {
        DBUG("Create GroundwaterFlowProcess.");

        // Process variable.
        {
            // Find the corresponding process variable.
            std::string const name = config.get<std::string>("process_variable");

            auto variable = std::find_if(variables.cbegin(), variables.cend(),
                    [&name](ProcessVariable const& v) {
                        return v.getName() == name;
                    });

            if (variable == variables.end())
                ERR("Expected process variable \'%s\' not found in provided variables list.",
                    name.c_str());

            DBUG("Associate hydraulic_head with process variable \'%s\'.",
                name.c_str());
            _hydraulic_head = const_cast<ProcessVariable*>(&*variable);
        }

        // Hydraulic conductivity parameter.
        {
            // find hydraulic_conductivity in process config
            boost::optional<std::string> const name =
                config.get_optional<std::string>("hydraulic_conductivity");
            if (!name)
            {
                ERR("Could not find required tag hydraulic_conductivity in "
                    "the process config.");
                std::abort();
            }

            // find corresponding parameter by name
            auto const parameter =
                std::find_if(parameters.cbegin(), parameters.cend(),
                             [&name](std::unique_ptr<ParameterBase> const& p)
                             {
                                 return p->name == name;
                             });

            if (parameter == parameters.end())
            {
                ERR("Could not find required parameter config for \'%s\' "
                    "among read parameters.",
                    name->c_str());
                std::abort();
            }

            _hydraulic_conductivity =
                dynamic_cast<const Parameter<double, const MeshLib::Element&>*>(
                    parameter->get());
            if (!_hydraulic_conductivity)
            {
                ERR("The hydraulic conductivity parameter is of incompatible "
                    "type.");
                std::abort();
            }
        }
    }

    template <unsigned GlobalDim>
    void createLocalAssemblers()
    {
        DBUG("Create local assemblers.");
        // Populate the vector of local assemblers.
        _local_assemblers.resize(_mesh.getNElements());
        // Shape matrices initializer
        using LocalDataInitializer = AssemblerLib::LocalDataInitializer<
            GroundwaterFlow::LocalAssemblerDataInterface,
            GroundwaterFlow::LocalAssemblerData,
            typename GlobalSetup::MatrixType,
            typename GlobalSetup::VectorType,
            GlobalDim>;

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
                *_hydraulic_conductivity,
                _integration_order);

        DBUG("Create global assembler.");
        _global_assembler.reset(
            new GlobalAssembler(*_A, *_rhs, *_local_to_global_index_map));

        DBUG("Initialize boundary conditions.");
        MeshGeoToolsLib::MeshNodeSearcher& hydraulic_head_mesh_node_searcher =
            MeshGeoToolsLib::MeshNodeSearcher::getMeshNodeSearcher(
                _hydraulic_head->getMesh());

        _hydraulic_head->initializeDirichletBCs(
                hydraulic_head_mesh_node_searcher,
                *_local_to_global_index_map, 0,
                _dirichlet_bc.global_ids, _dirichlet_bc.values);

        //
        // Neumann boundary conditions.
        //
        {
            // Find mesh nodes.
            MeshGeoToolsLib::BoundaryElementsSearcher hydraulic_head_mesh_element_searcher(
                _hydraulic_head->getMesh(), hydraulic_head_mesh_node_searcher);

            // Create a neumann BC for the hydraulic head storing them in the
            // _neumann_bcs vector.
            _hydraulic_head->createNeumannBcs(
                    std::back_inserter(_neumann_bcs),
                    hydraulic_head_mesh_element_searcher,
                    _global_setup,
                    _integration_order,
                    *_local_to_global_index_map,
                    0,
                    *_mesh_subset_all_nodes);
        }

        for (auto bc : _neumann_bcs)
            bc->initialize(_global_setup, *_A, *_rhs, _mesh.getDimension());

    }

    void initialize() override
    {
        DBUG("Initialize GroundwaterFlowProcess.");

        DBUG("Construct dof mappings.");
        // Create single component dof in every of the mesh's nodes.
        _mesh_subset_all_nodes = new MeshLib::MeshSubset(_mesh, &_mesh.getNodes());

        // Collect the mesh subsets in a vector.
        _all_mesh_subsets.push_back(new MeshLib::MeshSubsets(_mesh_subset_all_nodes));

        _local_to_global_index_map.reset(
            new AssemblerLib::LocalToGlobalIndexMap(_all_mesh_subsets, AssemblerLib::ComponentOrder::BY_COMPONENT));

#ifdef USE_PETSC
        DBUG("Allocate global matrix, vectors, and linear solver.");
        MathLib::PETScMatrixOption mat_opt;
        const MeshLib::NodePartitionedMesh& pmesh =
            static_cast<const MeshLib::NodePartitionedMesh&>(_mesh);
        mat_opt.d_nz = pmesh.getMaximumNConnectedNodesToNode();
        mat_opt.o_nz = mat_opt.d_nz;
        const std::size_t num_unknowns =
            _local_to_global_index_map->dofSizeGlobal();
        _A.reset(_global_setup.createMatrix(num_unknowns, mat_opt));
#else
        DBUG("Compute sparsity pattern");
        _node_adjacency_table.createTable(_mesh.getNodes());

        DBUG("Allocate global matrix, vectors, and linear solver.");
        const std::size_t num_unknowns = _local_to_global_index_map->dofSize();
        _A.reset(_global_setup.createMatrix(num_unknowns));
#endif

        _x.reset(_global_setup.createVector(num_unknowns));
        _rhs.reset(_global_setup.createVector(num_unknowns));
        _linearSolver.reset(new typename GlobalSetup::LinearSolver(*_A, "gw_"));

        setInitialConditions(*_hydraulic_head);

        if (_mesh.getDimension()==1)
            createLocalAssemblers<1>();
        else if (_mesh.getDimension()==2)
            createLocalAssemblers<2>();
        else if (_mesh.getDimension()==3)
            createLocalAssemblers<3>();
        else
            assert(false);
    }

    void setInitialConditions(ProcessVariable const& variable)
    {
        std::size_t const n = _mesh.getNNodes();
        for (std::size_t i = 0; i < n; ++i)
        {
            MeshLib::Location const l(_mesh.getID(),
                                      MeshLib::MeshItemType::Node, i);
            auto const global_index =
                _local_to_global_index_map->getGlobalIndex(
                    l, 0);  // 0 is the component id.
            _x->set(global_index,
                   variable.getInitialConditionValue(*_mesh.getNode(i)));
        }
    }

    bool solve(const double /*delta_t*/) override
    {
        DBUG("Solve GroundwaterFlowProcess.");

        _A->setZero();
        MathLib::setMatrixSparsity(*_A, _node_adjacency_table);
        *_rhs = 0;   // This resets the whole vector.

        // Call global assembler for each local assembly item.
        _global_setup.execute(*_global_assembler, _local_assemblers);

#ifdef USE_PETSC
        std::vector<PetscInt> dbc_pos;
        std::vector<PetscScalar> dbc_value;
        auto const& mesh =
            static_cast<const MeshLib::NodePartitionedMesh&>(_mesh);

        for (std::size_t i = 0; i < _dirichlet_bc.global_ids.size(); i++)
        {
            if (mesh.isGhostNode(_dirichlet_bc.global_ids[i]))
                continue;

            dbc_pos.push_back(static_cast<PetscInt>(
                mesh.getGlobalNodeID(_dirichlet_bc.global_ids[i])));
            dbc_value.push_back(
                static_cast<PetscScalar>(_dirichlet_bc.values[i]));
        }
        MathLib::applyKnownSolution(*_A, *_rhs, *_x, dbc_pos, dbc_value);
#else

        // Call global assembler for each Neumann boundary local assembler.
        for (auto bc : _neumann_bcs)
            bc->integrate(_global_setup);

        // Apply known values from the Dirichlet boundary conditions.
        MathLib::applyKnownSolution(*_A, *_rhs, _dirichlet_bc.global_ids, _dirichlet_bc.values);

#endif
        _linearSolver->solve(*_rhs, *_x);

        return true;
    }

    void post(std::string const& file_name) override
    {
        DBUG("Postprocessing GroundwaterFlowProcess.");

#ifdef USE_PETSC // Note: this is only a test

        std::vector<PetscScalar>  u(_x->size());
        _x->getGlobalVector(&u[0]);

        int rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        if(rank == 0)
        {
            std::string output_file_name = file_name;
            boost::erase_last(output_file_name, ".vtu");
            output_file_name.append(".dat");
            std::ofstream os(output_file_name);
            os << "SCALARS HEAD double 1"<<endl;
            os << "LOOKUP_TABLE default"<< endl;
            for (std::size_t i = 0; i < u.size(); ++i)
                 os << u[i] << "\n";
            os.close();
        }
#else
        std::string const property_name = "Result";

        // Get or create a property vector for results.
        boost::optional<MeshLib::PropertyVector<double>&> result;
        if (_mesh.getProperties().hasPropertyVector(property_name))
        {
            result = _mesh.getProperties().template
                getPropertyVector<double>(property_name);
        }
        else
        {
            result = _mesh.getProperties().template
                createNewPropertyVector<double>(property_name,
                    MeshLib::MeshItemType::Node);
            result->resize(_x->size());
        }
        assert(result && result->size() == _x->size());

        // Copy result
        for (std::size_t i = 0; i < _x->size(); ++i)
            (*result)[i] = (*_x)[i];

        // Write output file
        FileIO::VtuInterface vtu_interface(&_mesh, vtkXMLWriter::Binary, true);
        vtu_interface.writeToFile(file_name);
#endif
    }

    void postTimestep(std::string const& file_name, const unsigned /*timestep*/) override
    {
        post(file_name);
    }

    ~GroundwaterFlowProcess()
    {
        for (auto p : _neumann_bcs)
            delete p;

        for (auto p : _local_assemblers)
            delete p;

        for (auto p : _all_mesh_subsets)
            delete p;

        delete _mesh_subset_all_nodes;
    }

private:
    ProcessVariable* _hydraulic_head = nullptr;

    Parameter<double, MeshLib::Element const&> const* _hydraulic_conductivity = nullptr;

    MeshLib::MeshSubset const* _mesh_subset_all_nodes = nullptr;
    std::vector<MeshLib::MeshSubsets*> _all_mesh_subsets;

    GlobalSetup _global_setup;
    std::unique_ptr<typename GlobalSetup::LinearSolver> _linearSolver;
    std::unique_ptr<typename GlobalSetup::MatrixType> _A;
    std::unique_ptr<typename GlobalSetup::VectorType> _rhs;
    std::unique_ptr<typename GlobalSetup::VectorType> _x;

    using LocalAssembler = GroundwaterFlow::LocalAssemblerDataInterface<
        typename GlobalSetup::MatrixType, typename GlobalSetup::VectorType>;

    std::vector<LocalAssembler*> _local_assemblers;

    using GlobalAssembler =
        AssemblerLib::VectorMatrixAssembler<
            typename GlobalSetup::MatrixType,
            typename GlobalSetup::VectorType>;


    std::unique_ptr<AssemblerLib::LocalToGlobalIndexMap> _local_to_global_index_map;

    std::unique_ptr<GlobalAssembler> _global_assembler;

    /// Global ids in the global matrix/vector where the dirichlet bc is
    /// imposed and their corresponding values.
    struct DirichletBC {
        std::vector<std::size_t> global_ids;
        std::vector<double> values;
    } _dirichlet_bc;

    std::vector<NeumannBc<GlobalSetup>*> _neumann_bcs;

    MeshLib::NodeAdjacencyTable _node_adjacency_table;
};

}   // namespace ProcessLib

#endif  // PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_
