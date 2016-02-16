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

#include <boost/optional.hpp>

#include "AssemblerLib/LocalAssemblerBuilder.h"
#include "AssemblerLib/LocalDataInitializer.h"

#include "GroundwaterFlowFEM.h"
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
class GroundwaterFlowProcess final
        : public Process<GlobalSetup>
{
public:
    using Matrix = typename GlobalSetup::MatrixType;
    using Vector = typename GlobalSetup::VectorType;


    GroundwaterFlowProcess(
        MeshLib::Mesh& mesh,
        ProcessVariable& variable,
        Parameter<double, MeshLib::Element const&> const&
            hydraulic_conductivity,
        boost::optional<BaseLib::ConfigTree>&& linear_solver_options)
        : Process<GlobalSetup>(mesh),
          _hydraulic_conductivity(hydraulic_conductivity)
    {
        this->_process_variables.emplace_back(variable);
        if (linear_solver_options)
            Process<GlobalSetup>::setLinearSolverOptions(
                std::move(*linear_solver_options));
    }

    template <unsigned GlobalDim>
    void createLocalAssemblers()
    {
        DBUG("Create local assemblers.");
        // Populate the vector of local assemblers.
        _local_assemblers.resize(this->_mesh.getNElements());
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
            initializer, *this->_local_to_global_index_map);

        DBUG("Calling local assembler builder for all mesh elements.");
        this->_global_setup.execute(
                local_asm_builder,
                this->_mesh.getElements(),
                _local_assemblers,
                _hydraulic_conductivity,
                this->_integration_order);
    }

    std::string getLinearSolverName() const override
    {
        return "gw_";
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

    bool assemble(const double /*delta_t*/) override
    {
        DBUG("Assemble GroundwaterFlowProcess.");

        *this->_rhs = 0;   // This resets the whole vector.

        // Call global assembler for each local assembly item.
        this->_global_setup.execute(*this->_global_assembler,
                                    _local_assemblers);

        return true;
    }

    ~GroundwaterFlowProcess()
    {
        for (auto p : _local_assemblers)
            delete p;
    }

    //! \name ODESystem interface
    //! @{

    void assemble(const double t, Vector const& x,
                  Matrix& M, Matrix& K, Vector& b) override
    {
        // TODO implement;

        DBUG("Assemble GroundwaterFlowProcess.");

        *this->_rhs = 0;   // This resets the whole vector.

        // Call global assembler for each local assembly item.
        this->_global_setup.execute(*this->_global_assembler,
                                    _local_assemblers);
    }

    bool isLinear() const override
    {
        return true;
    }

    std::size_t getNumEquations() const override
    {
        // TODO implement
        return 1;
    }

    //! @}

private:
    Parameter<double, MeshLib::Element const&> const& _hydraulic_conductivity;

    using LocalAssembler = GroundwaterFlow::LocalAssemblerDataInterface<
        typename GlobalSetup::MatrixType, typename GlobalSetup::VectorType>;

    std::vector<LocalAssembler*> _local_assemblers;
};

template <typename GlobalSetup>
std::unique_ptr<GroundwaterFlowProcess<GlobalSetup>>
createGroundwaterFlowProcess(
    MeshLib::Mesh& mesh,
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

    // Linear solver options
    auto linear_solver_options = config.getConfSubtreeOptional("linear_solver");

    return std::unique_ptr<GroundwaterFlowProcess<GlobalSetup>>{
        new GroundwaterFlowProcess<GlobalSetup>{mesh, process_variable,
                                                hydraulic_conductivity,
                                                std::move(linear_solver_options)}};
}
}   // namespace ProcessLib

#endif  // PROCESS_LIB_GROUNDWATERFLOWPROCESS_H_
