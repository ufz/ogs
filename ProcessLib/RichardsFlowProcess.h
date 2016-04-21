/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_RICHARDSFLOWPROCESS_H_
#define PROCESS_LIB_RICHARDSFLOWPROCESS_H_

#include <cassert>

#include <boost/optional.hpp>

#include "AssemblerLib/LocalAssemblerBuilder.h"
#include "AssemblerLib/LocalDataInitializer.h"

#include "RichardsFlowFEM.h"
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
class RichardsFlowProcess final
        : public Process<GlobalSetup>
{
    // TODO change "this->" to "Base::"
    using Base = Process<GlobalSetup>;

public:
    using GlobalMatrix = typename GlobalSetup::MatrixType;
    using GlobalVector = typename GlobalSetup::VectorType;


	RichardsFlowProcess(
        MeshLib::Mesh& mesh,
        typename Process<GlobalSetup>::NonlinearSolver& nonlinear_solver,
        std::unique_ptr<typename Process<GlobalSetup>::TimeDiscretization>&& time_discretization,
        ProcessVariable& variable,
        Parameter<double, MeshLib::Element const&> const&
            intrinsic_permeability,
		Parameter<double, MeshLib::Element const&> const&
			porosity,
		Parameter<double, MeshLib::Element const&> const&
		viscosity,
		bool const has_gravity,
		std::map<std::string,
		std::unique_ptr<MathLib::PiecewiseLinearInterpolation >> const&
		curves)
        : Process<GlobalSetup>(mesh, nonlinear_solver, std::move(time_discretization)),
		_intrinsic_permeability(intrinsic_permeability),
		_porosity(porosity),
		_viscosity(viscosity),
		_has_gravity(has_gravity),
		_curves(curves)
    {
        this->_process_variables.emplace_back(variable);

        if (dynamic_cast<NumLib::ForwardEuler<GlobalVector>*>(
                    &Base::getTimeDiscretization()) != nullptr)
        {
            ERR("RichardsFlowProcess can not be solved with the ForwardEuler"
                " time discretization scheme. Aborting");
            // Because the M matrix is not assembled. Thus, the linearized system
            // would be singular. The same applies to CrankNicolson with theta = 0.0,
            // but this case is not checked here.
            // Anyway, the RichardsFlowProcess shall be transferred to a simpler
            // ODESystemTag in the future.
            std::abort();
        }
    }

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
			RichardsFlow::LocalAssemblerDataInterface,
			RichardsFlow::LocalAssemblerData,
            typename GlobalSetup::MatrixType,
            typename GlobalSetup::VectorType,
            GlobalDim,
            Parameter<double, MeshLib::Element const&> const&,
			Parameter<double, MeshLib::Element const&> const&,
			Parameter<double, MeshLib::Element const&> const&,
			bool const&,
			std::map<std::string,
			std::unique_ptr<MathLib::PiecewiseLinearInterpolation >> const&>;

        LocalDataInitializer initializer;

        using LocalAssemblerBuilder =
            AssemblerLib::LocalAssemblerBuilder<
                MeshLib::Element,
                LocalDataInitializer>;

		LocalAssemblerBuilder local_asm_builder(initializer, dof_table);

        DBUG("Calling local assembler builder for all mesh elements.");
        this->_global_setup.transform(
                local_asm_builder,
                mesh.getElements(),
                _local_assemblers,
                integration_order,
                _intrinsic_permeability,
				_porosity,
				_viscosity,
				_has_gravity,
				_curves);
    }

	void createAssemblers(AssemblerLib::LocalToGlobalIndexMap const& dof_table,
		MeshLib::Mesh const& mesh,
		unsigned const integration_order) override 
    {
		DBUG("Create global assembler.");
		_global_assembler.reset(new GlobalAssembler(dof_table));

		auto const dim = mesh.getDimension();
		if (dim == 1)
			createLocalAssemblers<1>(dof_table, mesh, integration_order);
		else if (dim == 2)
			createLocalAssemblers<2>(dof_table, mesh, integration_order);
		else if (dim == 3)
			createLocalAssemblers<3>(dof_table, mesh, integration_order);
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
    Parameter<double, MeshLib::Element const&> const& _intrinsic_permeability;
	Parameter<double, MeshLib::Element const&> const& _porosity;
	Parameter<double, MeshLib::Element const&> const& _viscosity;
	std::map<std::string,
		std::unique_ptr<MathLib::PiecewiseLinearInterpolation >> const&
		_curves;
	bool const _has_gravity;

    using LocalAssembler = RichardsFlow::LocalAssemblerDataInterface<
        typename GlobalSetup::MatrixType, typename GlobalSetup::VectorType>;

	using GlobalAssembler = AssemblerLib::VectorMatrixAssembler<
		GlobalMatrix, GlobalVector, LocalAssembler,
		NumLib::ODESystemTag::FirstOrderImplicitQuasilinear>;

    void assembleConcreteProcess(const double t, GlobalVector const& x,
                                 GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) override
    {
        // TODO It looks like, with little work this entire method can be moved to the Process class.

        DBUG("Assemble RichardsFlowProcess.");

        // Call global assembler for each local assembly item.
		this->_global_setup.executeMemberDereferenced(
			*_global_assembler, &GlobalAssembler::assemble,
			_local_assemblers, t, x, M, K, b);
    }
	std::unique_ptr<GlobalAssembler> _global_assembler;
    std::vector<std::unique_ptr<LocalAssembler>> _local_assemblers;
};

template <typename GlobalSetup>
std::unique_ptr<RichardsFlowProcess<GlobalSetup>>
createRichardsFlowProcess(
    MeshLib::Mesh& mesh,
    typename Process<GlobalSetup>::NonlinearSolver& nonlinear_solver,
    std::unique_ptr<typename Process<GlobalSetup>::TimeDiscretization>&& time_discretization,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
	BaseLib::ConfigTree const& config, 
	std::map<std::string,
			std::unique_ptr<MathLib::PiecewiseLinearInterpolation >> const&
	curves)
{
    config.checkConfParam("type", "RICHARDS_FLOW");

    DBUG("Create RichardsFlowProcess.");

    // Process variable.
    ProcessVariable& process_variable =
        findProcessVariable(config, "process_variable", variables);
    DBUG("Associate water_pressure with process variable \'%s\'.",
         process_variable.getName().c_str());

    // Intrinsic permeability parameter.
    auto& intrinsic_permeability =
        findParameter<double, MeshLib::Element const&>(
            config, "intrinsic_permeability", parameters);

    DBUG("Use \'%s\' as intrinsic_permeability parameter.",
		intrinsic_permeability.name.c_str());

	// Porosity parameter.
	auto& porosity =
		findParameter<double, MeshLib::Element const&>(
			config, "porosity", parameters);

	DBUG("Use \'%s\' as porosity parameter.",
		porosity.name.c_str());

	auto& viscosity =
		findParameter<double, MeshLib::Element const&>(
			config, "viscosity", parameters);

	DBUG("Use \'%s\' as viscosity parameter.",
		porosity.name.c_str());

	auto grav = config.getConfParam<bool>("g");

    return std::unique_ptr<RichardsFlowProcess<GlobalSetup>>{
        new RichardsFlowProcess<GlobalSetup>{
            mesh, nonlinear_solver,std::move(time_discretization),
            process_variable,
			intrinsic_permeability,
			porosity,
			viscosity,
			grav,curves
    }};
}
}   // namespace ProcessLib

#endif  // PROCESS_LIB_RICHARDSFLOWPROCESS_H_
