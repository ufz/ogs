/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file
 *
 * Created on August 19, 2016, 1:38 PM
 */

#pragma once

#include <memory>

#include "LiquidFlowLocalAssembler.h"
#include "LiquidFlowData.h"
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
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        LiquidFlowData&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        std::unique_ptr<ProcessLib::SurfaceFluxData>&& surfaceflux);

    void computeSecondaryVariableConcrete(double const t, double const dt,
                                          GlobalVector const& x,
                                          GlobalVector const& x_dot,
                                          int const process_id) override;

    bool isLinear() const override { return true; }

    Eigen::Vector3d getFlux(std::size_t const element_id,
                            MathLib::Point3d const& p,
                            double const t,
                            std::vector<GlobalVector*> const& x) const override;

    void postTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                     const double t,
                                     const double dt,
                                     int const process_id) override;

private:
    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh, unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, double const dt,
                                 std::vector<GlobalVector*> const& x,
                                 std::vector<GlobalVector*> const& xdot,
                                 int const process_id, GlobalMatrix& M,
                                 GlobalMatrix& K, GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        GlobalVector const& xdot, const double dxdot_dx, const double dx_dx,
        int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
        GlobalMatrix& Jac) override;

    LiquidFlowData _process_data;

    std::vector<std::unique_ptr<LiquidFlowLocalAssemblerInterface>>
        _local_assemblers;

    std::unique_ptr<ProcessLib::SurfaceFluxData> _surfaceflux;
    MeshLib::PropertyVector<double>* _hydraulic_flow = nullptr;
};

}  // namespace LiquidFlow
}  // namespace ProcessLib
