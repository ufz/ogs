/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Core>

#include "MathLib/Point3d.h"
#include "NumLib/NumericsConfig.h"

namespace NumLib
{
class LocalToGlobalIndexMap;
}  // namespace NumLib

namespace ProcessLib
{

/*! Common interface for local assemblers
 * NumLib::ODESystemTag::FirstOrderImplicitQuasilinear ODE systems.
 *
 * \todo Generalize to other NumLib::ODESystemTag's.
 */
class LocalAssemblerInterface
{
public:
    virtual ~LocalAssemblerInterface() = default;

    void setInitialConditions(std::size_t const mesh_item_id,
                              NumLib::LocalToGlobalIndexMap const& dof_table,
                              GlobalVector const& x, double const t,
                              bool const use_monolithic_scheme,
                              int const process_id);

    virtual void initialize(std::size_t const mesh_item_id,
                            NumLib::LocalToGlobalIndexMap const& dof_table);

    virtual void preAssemble(double const /*t*/, double const /*dt*/,
                             std::vector<double> const& /*local_x*/){};

    virtual void assemble(double const t, double const dt,
                          std::vector<double> const& local_x,
                          std::vector<double> const& local_x_prev,
                          std::vector<double>& local_M_data,
                          std::vector<double>& local_K_data,
                          std::vector<double>& local_b_data);

    virtual void assembleForStaggeredScheme(double const t, double const dt,
                                            Eigen::VectorXd const& local_x,
                                            Eigen::VectorXd const& local_x_prev,
                                            int const process_id,
                                            std::vector<double>& local_M_data,
                                            std::vector<double>& local_K_data,
                                            std::vector<double>& local_b_data);

    virtual void assembleWithJacobian(double const t, double const dt,
                                      std::vector<double> const& local_x,
                                      std::vector<double> const& local_x_prev,
                                      std::vector<double>& local_M_data,
                                      std::vector<double>& local_K_data,
                                      std::vector<double>& local_b_data,
                                      std::vector<double>& local_Jac_data);

    virtual void assembleWithJacobianForStaggeredScheme(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_x_prev, int const process_id,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data);

    virtual void computeSecondaryVariable(
        std::size_t const mesh_item_id,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_tables,
        double const t, double const dt, std::vector<GlobalVector*> const& x,
        GlobalVector const& x_prev, int const process_id);

    virtual void preTimestep(std::size_t const mesh_item_id,
                             NumLib::LocalToGlobalIndexMap const& dof_table,
                             GlobalVector const& x, double const t,
                             double const delta_t);

    virtual void postTimestep(
        std::size_t const mesh_item_id,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_tables,
        std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev, double const t,
        double const dt, bool const use_monolithic_scheme,
        int const process_id);

    void postNonLinearSolver(std::size_t const mesh_item_id,
                             NumLib::LocalToGlobalIndexMap const& dof_table,
                             std::vector<GlobalVector*> const& x,
                             std::vector<GlobalVector*> const& x_prev,
                             double const t, double const dt,
                             bool const use_monolithic_scheme,
                             int const process_id);

    /// Computes the flux in the point \c p_local_coords that is given in local
    /// coordinates using the values from \c local_x.
    /// Fits to monolithic scheme.
    virtual Eigen::Vector3d getFlux(
        MathLib::Point3d const& /*p_local_coords*/,
        double const /*t*/,
        std::vector<double> const& /*local_x*/) const
    {
        return Eigen::Vector3d{};
    }

    /// Fits to staggered scheme.
    virtual Eigen::Vector3d getFlux(
        MathLib::Point3d const& /*p_local_coords*/,
        double const /*t*/,
        std::vector<std::vector<double>> const& /*local_xs*/) const
    {
        return Eigen::Vector3d{};
    }

private:
    virtual void setInitialConditionsConcrete(
        std::vector<double> const& /*local_x*/, double const /*t*/,
        bool const /*use_monolithic_scheme*/, int const /*process_id*/)
    {
    }

    virtual void initializeConcrete() {}

    virtual void preTimestepConcrete(std::vector<double> const& /*local_x*/,
                                     double const /*t*/, double const /*dt*/)
    {
    }

    virtual void postTimestepConcrete(Eigen::VectorXd const& /*local_x*/,
                                      Eigen::VectorXd const& /*local_x_prev*/,
                                      double const /*t*/, double const /*dt*/,
                                      bool const /*use_monolithic_scheme*/,
                                      int const /*process_id*/)
    {
    }

    virtual void postNonLinearSolverConcrete(
        std::vector<double> const& /*local_x*/,
        std::vector<double> const& /*local_x_prev*/, double const /*t*/,
        double const /*dt*/, bool const /*use_monolithic_scheme*/,
        int const /*process_id*/)
    {
    }

    virtual void computeSecondaryVariableConcrete(
        double const /*t*/,
        double const /*dt*/,
        Eigen::VectorXd const& /*local_x*/,
        Eigen::VectorXd const& /*local_x_prev*/)
    {
    }

    virtual void computeSecondaryVariableWithCoupledProcessConcrete(
        double const /*t*/, std::vector<std::vector<double>> const&
        /*coupled_local_solutions*/)
    {
    }
};

}  // namespace ProcessLib
