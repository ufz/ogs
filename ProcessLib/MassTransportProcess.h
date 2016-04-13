/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_MASSTRANSPORTPROCESS_H_
#define PROCESS_LIB_MASSTRANSPORTPROCESS_H_

#include <cassert>

#include <boost/optional.hpp>

#include "AssemblerLib/LocalAssemblerBuilder.h"
#include "AssemblerLib/LocalDataInitializer.h"

#include "MassTransportFEM.h"
#include "Process.h"

namespace MeshLib
{
    class Element;
    class Mesh;
    template <typename PROP_VAL_TYPE> class PropertyVector;
}

namespace ProcessLib
{

template<typename GlobalSetup>
class MassTransportProcess final
        : public Process<GlobalSetup>
{
    // TODO change "this->" to "Base::"
    using Base = Process<GlobalSetup>;

public:
    using GlobalMatrix = typename GlobalSetup::MatrixType;
    using GlobalVector = typename GlobalSetup::VectorType;

	//template <unsigned GlobalDim>
    MassTransportProcess(
        MeshLib::Mesh& mesh,
        typename Process<GlobalSetup>::NonlinearSolver& nonlinear_solver,
        std::unique_ptr<typename Process<GlobalSetup>::TimeDiscretization>&& time_discretization,
        ProcessVariable& variable,
        Parameter<double, MeshLib::Element const&> const&
            diffusion_coefficient,
		Parameter<double, MeshLib::Element const&> const&
			velocity)
        : Process<GlobalSetup>(mesh, nonlinear_solver, std::move(time_discretization)),
          _diffusion_coefficient(diffusion_coefficient),
		_velocity(velocity)
    {
        this->_process_variables.emplace_back(variable);

        if (dynamic_cast<NumLib::ForwardEuler<GlobalVector>*>(
                    &Base::getTimeDiscretization()) != nullptr)
        {
            ERR("MassTransportProcess can not be solved with the ForwardEuler"
                " time discretization scheme. Aborting");
            // Because the M matrix is not assembled. Thus, the linearized system
            // would be singular. The same applies to CrankNicolson with theta = 0.0,
            // but this case is not checked here.
            // Anyway, the GroundwaterFlowProcess shall be transferred to a simpler
            // ODESystemTag in the future.
            std::abort();
        }
    }

    template <unsigned GlobalDim>
    void createLocalAssemblers()
    {
        DBUG("Create local assemblers.");
        // Populate the vector of local assemblers.
        _local_assemblers.resize(this->_mesh.getNElements());
        // Shape matrices initializer
        using LocalDataInitializer = AssemblerLib::LocalDataInitializer<
            MassTransport::LocalAssemblerDataInterface,
            MassTransport::LocalAssemblerData,
            typename GlobalSetup::MatrixType,
            typename GlobalSetup::VectorType,
            GlobalDim,
            Parameter<double, MeshLib::Element const&> const&,
			Parameter<double, MeshLib::Element const&> const&>;

        LocalDataInitializer initializer;

        using LocalAssemblerBuilder =
            AssemblerLib::LocalAssemblerBuilder<
                MeshLib::Element,
                LocalDataInitializer>;

        LocalAssemblerBuilder local_asm_builder(
            initializer, *this->_local_to_global_index_map);

        DBUG("Calling local assembler builder for all mesh elements.");
        this->_global_setup.transform(
                local_asm_builder,
                this->_mesh.getElements(),
                _local_assemblers,
                this->_integration_order,
                _diffusion_coefficient,
			_velocity);
    }

    void createLocalAssemblers() override
    {
        if (this->_mesh.getDimension()==1)
            createLocalAssemblers<1>();
        else if (this->_mesh.getDimension()==2)
            createLocalAssemblers<2>();
        else if (this->_mesh.getDimension()==3)
            createLocalAssemblers<3>();
        else
            assert(false);
    }

    //! \name ODESystem interface
    //! @{

    bool isLinear() const override
    {
        return true;
    }

    //! @}

private:
    Parameter<double, MeshLib::Element const&> const& _diffusion_coefficient;
	Parameter<double, MeshLib::Element const&> const& _velocity;
    using LocalAssembler = MassTransport::LocalAssemblerDataInterface<
        typename GlobalSetup::MatrixType, typename GlobalSetup::VectorType>;


    void assembleConcreteProcess(const double t, GlobalVector const& x,
                                 GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) override
    {
        // TODO It looks like, with little work this entire method can be moved to the Process class.

        DBUG("Assemble GroundwaterFlowProcess.");

        // Call global assembler for each local assembly item.
        this->_global_setup.executeDereferenced(
            *this->_global_assembler, _local_assemblers, t, x, M, K, b);
    }

    std::vector<std::unique_ptr<LocalAssembler>> _local_assemblers;
};

template <typename GlobalSetup>
std::unique_ptr<MassTransportProcess<GlobalSetup>>
createMassTransportProcess(
    MeshLib::Mesh& mesh,
    typename Process<GlobalSetup>::NonlinearSolver& nonlinear_solver,
    std::unique_ptr<typename Process<GlobalSetup>::TimeDiscretization>&& time_discretization,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    config.checkConfParam("type", "Mass_Transport");

    DBUG("Create MassTransportProcess.");

    // Process variable.
    ProcessVariable& process_variable =
        findProcessVariable(config, "process_variable", variables);
    DBUG("Associate concentration with process variable \'%s\'.",
         process_variable.getName().c_str());

    // Diffusion coefficient parameter.
    auto& diffusion_coefficient =
        findParameter<double, MeshLib::Element const&>(
            config, "diffusion_coefficient", parameters);
	DBUG("Use \'%s\' as diffusion coefficient parameter.",
         diffusion_coefficient.name.c_str());
	
	//velocity parameter
	auto& velocity =
		findParameter<double, MeshLib::Element const&>(
			config, "velocity", parameters);
	DBUG("Use \'%s\' as velocity parameter.",
		velocity.name.c_str());

    return std::unique_ptr<MassTransportProcess<GlobalSetup>>{
        new MassTransportProcess<GlobalSetup>{
            mesh, nonlinear_solver,std::move(time_discretization),
            process_variable,
			diffusion_coefficient,
			velocity
    }};
}
}   // namespace ProcessLib

#endif  // PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_
