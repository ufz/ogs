/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   LiquidFlowProcess.h
 *
 * Created on August 19, 2016, 1:38 PM
 */

#pragma once

#include <memory>

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/Process.h"

#include "MaterialLib/Fluid/FluidProperties/FluidProperties.h"

namespace MeshLib
{
class Element;
class Mesh;
template <typename PROP_VAL_TYPE>
class PropertyVector;
}

namespace ProcessLib
{
namespace LiquidFlow
{
class LiquidFlowLocalAssemblerInterface;
class LiquidFlowMaterialProperties;

/**
 * \brief A class to simulate the liquid flow process in porous media described
 * by
 *
 * \f[
 *     \frac{\partial n \rho_l}{\partial T} \frac{\partial T}{\partial t}/\rho_l
 *       + (\frac{\partial n \rho_l}{\partial p}/\rho_l + \beta_s)
 *          \frac{\partial p}{\partial t}
 *       -\nabla (\frac{K}{\mu}(\nabla p + \rho_l g \nabla z) ) = Q
 * \f]
 * where
 *    \f{eqnarray*}{
 *       &p:&        \mbox{pore pressure,}\\
 *       &T: &       \mbox{Temperature,}\\
 *       &\rho_l:&   \mbox{liquid density,}\\
 *       &\beta_s:&  \mbox{specific storage,}\\
 *       &K:&        \mbox{permeability,}\\
 *       &\mu:&      \mbox{viscosity,}\\
 *    \f}
 */
class LiquidFlowProcess final : public Process
{
public:
    LiquidFlowProcess(
        MeshLib::Mesh& mesh,
        std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        unsigned const integration_order,
        std::vector<std::reference_wrapper<ProcessVariable>>&&
            process_variables,
        SecondaryVariableCollection&& secondary_variables,
        NumLib::NamedFunctionCaller&& named_function_caller,
        MeshLib::PropertyVector<int> const& material_ids,
        bool const has_material_ids,
        int const gravitational_axis_id,
        double const gravitational_acceleration,
        double const reference_temperature,
        BaseLib::ConfigTree const& config);

    void computeSecondaryVariableConcrete(
        double const t,
        GlobalVector const& x,
        StaggeredCouplingTerm const& coupled_term) override;

    bool isLinear() const override { return true; }
    int getGravitationalAxisID() const { return _gravitational_axis_id; }
    double getGravitationalAcceleration() const
    {
        return _gravitational_acceleration;
    }

    LiquidFlowMaterialProperties* getLiquidFlowMaterialProperties() const
    {
        return _material_properties.get();
    }

private:
    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh, unsigned const integration_order) override;

    void assembleConcreteProcess(
        const double t, GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
        GlobalVector& b, StaggeredCouplingTerm const& coupling_term) override;

    void assembleWithJacobianConcreteProcess(
        const double t, GlobalVector const& x, GlobalVector const& xdot,
        const double dxdot_dx, const double dx_dx, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac,
        StaggeredCouplingTerm const& coupling_term) override;

    const int _gravitational_axis_id;
    const double _gravitational_acceleration;
    const double _reference_temperature;
    const std::unique_ptr<LiquidFlowMaterialProperties> _material_properties;

    std::vector<std::unique_ptr<LiquidFlowLocalAssemblerInterface>>
        _local_assemblers;
};

}  // end of namespace
}  // end of namespace
