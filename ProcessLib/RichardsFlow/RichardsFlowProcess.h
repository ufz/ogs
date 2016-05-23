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

#include "AssemblerLib/VectorMatrixAssembler.h"
#include "NumLib/Extrapolation/LocalLinearLeastSquaresExtrapolator.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"

#include "RichardsFlowFEM.h"
#include "RichardsFlowProcessData.h"

namespace MeshLib
{
    class Element;
}

namespace ProcessLib
{

namespace RichardsFlow
{

template<typename GlobalSetup>
class RichardsFlowProcess final
        : public Process<GlobalSetup>
{
    using Base = Process<GlobalSetup>;
    using GlobalMatrix = typename GlobalSetup::MatrixType;
    using GlobalVector = typename GlobalSetup::VectorType;

public:
	RichardsFlowProcess(
        MeshLib::Mesh& mesh,
        typename Base::NonlinearSolver& nonlinear_solver,
        std::unique_ptr<typename Base::TimeDiscretization>&& time_discretization,
        std::vector<std::reference_wrapper<ProcessVariable>>&& process_variables,
		RichardsFlowProcessData&& process_data,
		SecondaryVariableCollection<GlobalVector>&& secondary_variables,
		ProcessOutput<GlobalVector>&& process_output
		)
        : Process<GlobalSetup>(mesh, nonlinear_solver, std::move(time_discretization),
                               std::move(process_variables),
							   std::move(secondary_variables),
							   std::move(process_output))
        , _process_data(std::move(process_data))
    {
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

    //! \name ODESystem interface
    //! @{

    bool isLinear() const override
    {
        return false;
    }

    //! @}

private:
    using LocalAssemblerInterface =
        RichardsFlowLocalAssemblerInterface<GlobalMatrix, GlobalVector>;

    using GlobalAssembler = AssemblerLib::VectorMatrixAssembler<
            GlobalMatrix, GlobalVector, LocalAssemblerInterface,
            NumLib::ODESystemTag::FirstOrderImplicitQuasilinear>;

	using ExtrapolatorInterface = NumLib::Extrapolator<
		GlobalVector, IntegrationPointValue, LocalAssemblerInterface>;
	using ExtrapolatorImplementation = NumLib::LocalLinearLeastSquaresExtrapolator<
		GlobalVector, IntegrationPointValue, LocalAssemblerInterface>;

	virtual void preIteration(const unsigned iter, GlobalVector const& x) override
	{
		//(void)iter; (void)x; // by default do nothing
		//x.write("test" + std::to_string(iter) );//
	}
	void initializeConcreteProcess(
		AssemblerLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override
    {
        DBUG("Create global assembler.");
        _global_assembler.reset(new GlobalAssembler(dof_table));

        ProcessLib::createLocalAssemblers<GlobalSetup, LocalAssemblerData>(
                    mesh.getDimension(), mesh.getElements(),
                    dof_table, integration_order, _local_assemblers,
                    _process_data);

		_extrapolator.reset(
			new ExtrapolatorImplementation(Base::getMatrixSpecifications()));

		Base::_secondary_variables.addSecondaryVariable(
			"saturation", 1,
			makeExtrapolator(IntegrationPointValue::Saturation, *_extrapolator,
			_local_assemblers));
    }

    void assembleConcreteProcess(const double t, GlobalVector const& x,
                                 GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) override
    {
        DBUG("Assemble RichardsFlowProcess.");

        // Call global assembler for each local assembly item.
        GlobalSetup::executeMemberDereferenced(
            *_global_assembler, &GlobalAssembler::assemble,
            _local_assemblers, t, x, M, K, b);
    }


	RichardsFlowProcessData _process_data;

    std::unique_ptr<GlobalAssembler> _global_assembler;
    std::vector<std::unique_ptr<LocalAssemblerInterface>> _local_assemblers;

	std::unique_ptr<ExtrapolatorInterface> _extrapolator;
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
    auto process_variables =
        findProcessVariables(variables, config, { "process_variable" });

	//

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

	// Viscosity parameter.
	auto& viscosity =
		findParameter<double, MeshLib::Element const&>(
			config, "viscosity", parameters);

	DBUG("Use \'%s\' as viscosity parameter.",
		viscosity.name.c_str());

	//has  gravity
	auto grav = config.getConfParam<bool>("g");

    RichardsFlowProcessData process_data {
		intrinsic_permeability,
		porosity,
		viscosity,
		grav,
		curves
    };

	SecondaryVariableCollection<typename GlobalSetup::VectorType>
		secondary_variables{ config.getConfSubtreeOptional("secondary_variables"),
		{ "saturation"} };

	ProcessOutput<typename GlobalSetup::VectorType>
		process_output{ config.getConfSubtree("output"),
		process_variables, secondary_variables };

    return std::unique_ptr<RichardsFlowProcess<GlobalSetup>>{
        new RichardsFlowProcess<GlobalSetup>{
            mesh, nonlinear_solver,std::move(time_discretization),
            std::move(process_variables),
            std::move(process_data),
			std::move(secondary_variables),
			std::move(process_output)
        }
    };
}

}   // namespace RichardsFlow
}   // namespace ProcessLib

#endif  // PROCESS_LIB_RICHARDSFLOWPROCESS_H_
