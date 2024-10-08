/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * Created on August 19, 2016, 1:38 PM
 */

#pragma once

#include <memory>

#include "LiquidFlowData.h"
#include "LiquidFlowLocalAssembler.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/SurfaceFlux/SurfaceFluxData.h"

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
/**
 * \brief A class to simulate the single phase flow process
 *  in porous media described by the governing equation:
 *
 * \f[
 *        \left(\frac{\partial}{\partial p}(\phi \rho) + \rho \beta_s\right)
 *          \frac{\partial p}{\partial t}
 *       -\nabla \left(\rho\frac{\mathbf K}{\mu}(\nabla p + \rho g \nabla
 * z)\right ) = \rho Q, \f] where \f{eqnarray*}{
 *       &p:&       \mbox{pore pressure,}\\
 *       &\phi:&    \mbox{porosity,}\\
 *       &\rho:&    \mbox{liquid or gas density,}\\
 *       &\beta_s:& \mbox{specific storage,}\\
 *       &{\mathbf K}:&  \mbox{permeability,}\\
 *       &\mu:&     \mbox{viscosity,}\\
 *       &g:&       \mbox{gravitational constant,}\\
 *       &Q:&       \mbox{Source/sink term in m}^3/{\text s}.\\
 *    \f}
 * This governing equation represents the mass balance.
 *
 * If the density is assumed constant, for example for a groundwater modelling,
 *  the governing equation is scaled with the density, and it becomes volume
 *  balanced as:
 * \f[
 *        \left(\frac{1}{\rho}\frac{\partial}{\partial p}(\phi \rho)
 *          +  \beta_s\right)
 *          \frac{\partial p}{\partial t}
 *       -\nabla \left(\frac{\mathbf K}{\mu}(\nabla p + \rho g \nabla z)\right )
 *      =  Q,
 * \f]
 *
 * An optional input tag `equation_balance_type` of this process can be used to
 * select whether to use the volume balanced equation or the mass balanced
 * equation. <b>By default, we assume that volume balanced equation is used</b>.
 *
 * Be aware that the Neumann condition is
 *    \f{eqnarray*}{
 *       & -\frac{\mathbf K}{\mu}(\nabla p + \rho g \nabla z)
 *          \cdot \mathbf n = q_v [\text{m/s}]: &
 *          \mbox{ for the volume balance equation,}\\
 *       & -\rho\frac{\mathbf K}{\mu}(\nabla p + \rho g \nabla z)
 *          \cdot \mathbf n = q_f [\text{kg/m²/s}]: &
 *            \mbox{for the mass balance equation,}
 *    \f}
 *   with \f$ \mathbf n \f$ the outer normal of the boundary.
 *  */
class LiquidFlowProcess final : public Process
{
public:
    LiquidFlowProcess(
        std::string name, MeshLib::Mesh& mesh,
        std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        LiquidFlowData&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        std::unique_ptr<ProcessLib::SurfaceFluxData>&& surfaceflux,
        bool const is_linear);

    void computeSecondaryVariableConcrete(double const t, double const dt,
                                          std::vector<GlobalVector*> const& x,
                                          GlobalVector const& x_prev,
                                          int const process_id) override;

    bool isLinear() const override { return _is_linear; }

    Eigen::Vector3d getFlux(std::size_t const element_id,
                            MathLib::Point3d const& p,
                            double const t,
                            std::vector<GlobalVector*> const& x) const override;

    void postTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                     std::vector<GlobalVector*> const& x_prev,
                                     const double t,
                                     const double dt,
                                     int const process_id) override;

private:
    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh, unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, double const dt,
                                 std::vector<GlobalVector*> const& x,
                                 std::vector<GlobalVector*> const& x_prev,
                                 int const process_id, GlobalMatrix& M,
                                 GlobalMatrix& K, GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev, int const process_id,
        GlobalVector& b, GlobalMatrix& Jac) override;

    LiquidFlowData _process_data;

    std::vector<std::unique_ptr<LiquidFlowLocalAssemblerInterface>>
        _local_assemblers;

    std::unique_ptr<ProcessLib::SurfaceFluxData> _surfaceflux;
    MeshLib::PropertyVector<double>* _hydraulic_flow = nullptr;
    bool _is_linear = false;
};

}  // namespace LiquidFlow
}  // namespace ProcessLib
