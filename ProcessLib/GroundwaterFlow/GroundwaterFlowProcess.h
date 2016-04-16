/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_
#define PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_

#include <cassert>

#include "AssemblerLib/LocalAssemblerBuilder.h"
#include "AssemblerLib/LocalDataInitializer.h"
#include "ProcessLib/Process.h"

#include "GroundwaterFlowFEM.h"

namespace MeshLib
{
    class Element;
    class Mesh;
    template <typename PROP_VAL_TYPE> class PropertyVector;
}

namespace ProcessLib
{

namespace GroundwaterFlow
{

template<typename GlobalSetup>
class GroundwaterFlowProcess final
        : public Process<GlobalSetup>
{
    using Base = Process<GlobalSetup>;

public:
    using GlobalMatrix = typename GlobalSetup::MatrixType;
    using GlobalVector = typename GlobalSetup::VectorType;


    GroundwaterFlowProcess(
        MeshLib::Mesh& mesh,
        typename Process<GlobalSetup>::NonlinearSolver& nonlinear_solver,
        std::unique_ptr<typename Process<GlobalSetup>::TimeDiscretization>&& time_discretization,
        ProcessVariable& variable,
        Parameter<double, MeshLib::Element const&> const&
            hydraulic_conductivity)
        : Process<GlobalSetup>(mesh, nonlinear_solver, std::move(time_discretization)),
          _hydraulic_conductivity(hydraulic_conductivity)
    {
        Base::_process_variables.emplace_back(variable);

        if (dynamic_cast<NumLib::ForwardEuler<GlobalVector>*>(
                    &Base::getTimeDiscretization()) != nullptr)
        {
            ERR("GroundwaterFlowProcess can not be solved with the ForwardEuler"
                " time discretization scheme. Aborting");
            // Because the M matrix is not assembled. Thus, the linearized system
            // would be singular. The same applies to CrankNicolson with theta = 0.0,
            // but this case is not checked here.
            // Anyway, the GroundwaterFlowProcess shall be transferred to a simpler
            // ODESystemTag in the future.
            std::abort();
        }
    }

    //! \name ODESystem interface
    //! @{

    bool isLinear() const override
    {
        return true;
    }

    //! @}

private:
    Parameter<double, MeshLib::Element const&> const& _hydraulic_conductivity;

    using LocalAssemblerInterface =
        ProcessLib::LocalAssemblerInterface<GlobalMatrix, GlobalVector>;

    using GlobalAssembler = AssemblerLib::VectorMatrixAssembler<
            GlobalMatrix, GlobalVector, LocalAssemblerInterface,
            NumLib::ODESystemTag::FirstOrderImplicitQuasilinear>;

    template <unsigned GlobalDim>
    void createLocalAssemblers(AssemblerLib::LocalToGlobalIndexMap const& dof_table,
                               MeshLib::Mesh const& mesh,
                               unsigned const integration_order)
    {
        DBUG("Create local assemblers.");
        // Populate the vector of local assemblers.
        _local_assemblers.resize(mesh.getNElements());
        // Shape matrices initializer
        using LocalDataInitializer = AssemblerLib::LocalDataInitializer<
            ProcessLib::LocalAssemblerInterface,
            GroundwaterFlow::LocalAssemblerData,
            typename GlobalSetup::MatrixType,
            typename GlobalSetup::VectorType,
            GlobalDim,
            ProcessLib::Parameter<double, MeshLib::Element const&> const&>;

        LocalDataInitializer initializer;

        using LocalAssemblerBuilder =
            AssemblerLib::LocalAssemblerBuilder<
                MeshLib::Element,
                LocalDataInitializer>;

        LocalAssemblerBuilder local_asm_builder(initializer, dof_table);

        DBUG("Calling local assembler builder for all mesh elements.");
        GlobalSetup::transform(
                local_asm_builder,
                mesh.getElements(),
                _local_assemblers,
                integration_order,
                _hydraulic_conductivity);
    }

    void createAssemblers(AssemblerLib::LocalToGlobalIndexMap const& dof_table,
                          MeshLib::Mesh const& mesh,
                          unsigned const integration_order) override
    {
        DBUG("Create global assembler.");
        _global_assembler.reset(new GlobalAssembler(dof_table));

        auto const dim = mesh.getDimension();
        if (dim==1)
            createLocalAssemblers<1>(dof_table, mesh, integration_order);
        else if (dim==2)
            createLocalAssemblers<2>(dof_table, mesh, integration_order);
        else if (dim==3)
            createLocalAssemblers<3>(dof_table, mesh, integration_order);
        else
            assert(false);
    }

    void assembleConcreteProcess(const double t, GlobalVector const& x,
                                 GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) override
    {
        // TODO It looks like, with little work this entire method can be moved to the Process class.

        DBUG("Assemble GroundwaterFlowProcess.");

        // Call global assembler for each local assembly item.
        GlobalSetup::executeMemberDereferenced(
            *_global_assembler, &GlobalAssembler::assemble,
            _local_assemblers, t, x, M, K, b);
    }

    std::unique_ptr<GlobalAssembler> _global_assembler;
    std::vector<std::unique_ptr<LocalAssemblerInterface>> _local_assemblers;
};

template <typename GlobalSetup>
std::unique_ptr<GroundwaterFlowProcess<GlobalSetup>>
createGroundwaterFlowProcess(
    MeshLib::Mesh& mesh,
    typename Process<GlobalSetup>::NonlinearSolver& nonlinear_solver,
    std::unique_ptr<typename Process<GlobalSetup>::TimeDiscretization>&& time_discretization,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    config.checkConfParam("type", "GROUNDWATER_FLOW");

    DBUG("Create GroundwaterFlowProcess.");

    // Process variable.
    ProcessVariable& process_variable =
        findProcessVariable(config, "process_variable", variables);
    DBUG("Associate hydraulic_head with process variable \'%s\'.",
         process_variable.getName().c_str());

    // Hydraulic conductivity parameter.
    auto& hydraulic_conductivity =
        findParameter<double, MeshLib::Element const&>(
            config, "hydraulic_conductivity", parameters);

    DBUG("Use \'%s\' as hydraulic conductivity parameter.",
         hydraulic_conductivity.name.c_str());

    return std::unique_ptr<GroundwaterFlowProcess<GlobalSetup>>{
        new GroundwaterFlowProcess<GlobalSetup>{
            mesh, nonlinear_solver,std::move(time_discretization),
            process_variable,
            hydraulic_conductivity
    }};
}

}   // namespace GroundwaterFlow
}   // namespace ProcessLib

#endif  // PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_
