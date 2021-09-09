/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <vector>

#include "IntegrationPointData.h"
#include "LocalAssemblerInterface.h"
#include "MaterialLib/MPL/VariableType.h"
#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MathLib/KelvinVector.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "RichardsMechanicsProcessData.h"

namespace ProcessLib
{
namespace RichardsMechanics
{
namespace MPL = MaterialPropertyLib;

/// Used for the extrapolation of the integration point values. It is ordered
/// (and stored) by integration points.
template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType, Eigen::aligned_allocator<ShapeMatrixType>> N_u;
};

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
class RichardsMechanicsLocalAssembler
    : public LocalAssemblerInterface<DisplacementDim>
{
public:
    using ShapeMatricesTypeDisplacement =
        ShapeMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;
    using ShapeMatricesTypePressure =
        ShapeMatrixPolicyType<ShapeFunctionPressure, DisplacementDim>;

    using GlobalDimMatrixType =
        typename ShapeMatricesTypePressure::GlobalDimMatrixType;

    using BMatricesType =
        BMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;
    using KelvinVectorType = typename BMatricesType::KelvinVectorType;

    using IpData =
        IntegrationPointData<BMatricesType, ShapeMatricesTypeDisplacement,
                             ShapeMatricesTypePressure, DisplacementDim,
                             ShapeFunctionDisplacement::NPOINTS>;

    static int const KelvinVectorSize =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;

    using SymmetricTensor = Eigen::Matrix<double, KelvinVectorSize, 1>;

    RichardsMechanicsLocalAssembler(RichardsMechanicsLocalAssembler const&) =
        delete;
    RichardsMechanicsLocalAssembler(RichardsMechanicsLocalAssembler&&) = delete;

    RichardsMechanicsLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        RichardsMechanicsProcessData<DisplacementDim>& process_data);

    /// \return the number of read integration points.
    std::size_t setIPDataInitialConditions(
        std::string const& name,
        double const* values,
        int const integration_order) override;

    void setInitialConditionsConcrete(std::vector<double> const& local_x,
                                      double const t,
                                      bool const use_monolithic_scheme,
                                      int const process_id) override;

    void assemble(double const t, double const dt,
                  std::vector<double> const& local_x,
                  std::vector<double> const& local_xdot,
                  std::vector<double>& local_M_data,
                  std::vector<double>& local_K_data,
                  std::vector<double>& local_rhs_data) override;

    void assembleWithJacobian(double const t, double const dt,
                              std::vector<double> const& local_x,
                              std::vector<double> const& local_xdot,
                              const double /*dxdot_dx*/, const double /*dx_dx*/,
                              std::vector<double>& /*local_M_data*/,
                              std::vector<double>& /*local_K_data*/,
                              std::vector<double>& local_rhs_data,
                              std::vector<double>& local_Jac_data) override;

    void assembleWithJacobianForStaggeredScheme(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_xdot, const double dxdot_dx,
        const double dx_dx, int const process_id,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data,
        std::vector<double>& local_Jac_data) override;

    void initializeConcrete() override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto& ip_data = _ip_data[ip];

            /// Set initial stress from parameter.
            if (_process_data.initial_stress != nullptr)
            {
                ParameterLib::SpatialPosition const x_position{
                    std::nullopt, _element.getID(), ip,
                    MathLib::Point3d(NumLib::interpolateCoordinates<
                                     ShapeFunctionDisplacement,
                                     ShapeMatricesTypeDisplacement>(
                        _element, ip_data.N_u))};

                ip_data.sigma_eff =
                    MathLib::KelvinVector::symmetricTensorToKelvinVector<
                        DisplacementDim>((*_process_data.initial_stress)(
                        std::numeric_limits<
                            double>::quiet_NaN() /* time independent */,
                        x_position));
            }

            ip_data.pushBackState();
        }
    }

    void postTimestepConcrete(Eigen::VectorXd const& /*local_x*/,
                              double const /*t*/,
                              double const /*dt*/) override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data[ip].pushBackState();
        }
    }

    void computeSecondaryVariableConcrete(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_x_dot) override;

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N_u = _secondary_data.N_u[integration_point];

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N_u.data(), N_u.size());
    }

    std::vector<double> getSigma() const override;

    std::vector<double> const& getIntPtDarcyVelocity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override;

    std::vector<double> getSaturation() const override;
    std::vector<double> const& getIntPtSaturation(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override;

    std::vector<double> getMicroSaturation() const override;
    std::vector<double> const& getIntPtMicroSaturation(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override;

    std::vector<double> getMicroPressure() const override;
    std::vector<double> const& getIntPtMicroPressure(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override;

    std::vector<double> getPorosity() const override;
    std::vector<double> const& getIntPtPorosity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override;

    std::vector<double> getTransportPorosity() const override;
    std::vector<double> const& getIntPtTransportPorosity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override;

    std::vector<double> const& getIntPtSigma(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override;

    std::vector<double> getSwellingStress() const override;
    std::vector<double> const& getIntPtSwellingStress(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override;

    std::vector<double> getEpsilon() const override;
    std::vector<double> const& getIntPtEpsilon(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override;

    std::vector<double> getMaterialStateVariableInternalState(
        std::function<BaseLib::DynamicSpan<double>(
            typename MaterialLib::Solids::MechanicsBase<DisplacementDim>::
                MaterialStateVariables&)> const& get_values_span,
        int const& n_components) const override;

    std::vector<double> const& getIntPtDryDensitySolid(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override;

private:
    /**
     * Assemble local matrices and vectors arise from the linearized discretized
     * weak form of the residual of the momentum balance equation,
     *      \f[
     *            \nabla (\sigma - \alpha_b p \mathrm{I}) = f
     *      \f]
     * where \f$ \sigma\f$ is the effective stress tensor, \f$p\f$ is the pore
     * pressure, \f$\alpha_b\f$ is the Biot constant, \f$\mathrm{I}\f$ is the
     * identity tensor, and \f$f\f$ is the body force.
     *
     * @param t               Time
     * @param dt              Time increment
     * @param local_x         Nodal values of \f$x\f$ of an element.
     * @param local_xdot      Nodal values of \f$\dot{x}\f$ of an element.
     * @param dxdot_dx        Value of \f$\dot{x} \cdot dx\f$.
     * @param dx_dx           Value of \f$ x \cdot dx\f$.
     * @param local_M_data    Mass matrix of an element, which takes the form of
     *                        \f$ \int N^T N\mathrm{d}\Omega\f$. Not used.
     * @param local_K_data    Laplacian matrix of an element, which takes the
     *         form of \f$ \int (\nabla N)^T K \nabla N\mathrm{d}\Omega\f$.
     *                        Not used.
     * @param local_b_data    Right hand side vector of an element.
     * @param local_Jac_data  Element Jacobian matrix for the Newton-Raphson
     *                        method.
     */
    void assembleWithJacobianForDeformationEquations(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_xdot, const double dxdot_dx,
        const double dx_dx, std::vector<double>& local_M_data,
        std::vector<double>& local_K_data, std::vector<double>& local_b_data,
        std::vector<double>& local_Jac_data);

    /**
     * Assemble local matrices and vectors arise from the linearized discretized
     * weak form of the residual of the mass balance equation of single phase
     * flow,
     *      \f[
     *          \alpha \cdot{p} - \nabla (K (\nabla p + \rho g \nabla z) +
     *          \alpha_b \nabla \cdot \dot{u}  = Q
     *      \f]
     * where \f$ alpha\f$ is a coefficient may depend on storage or the fluid
     * density change, \f$ \rho\f$ is the fluid density, \f$g\f$ is the
     * gravitational acceleration, \f$z\f$ is the vertical coordinate, \f$u\f$
     * is the displacement, and \f$Q\f$ is the source/sink term.
     *
     * @param t               Time
     * @param dt              Time increment
     * @param local_x         Nodal values of \f$x\f$ of an element.
     * @param local_xdot      Nodal values of \f$\dot{x}\f$ of an element.
     * @param dxdot_dx        Value of \f$\dot{x} \cdot dx\f$.
     * @param dx_dx           Value of \f$ x \cdot dx\f$.
     * @param local_M_data    Mass matrix of an element, which takes the form of
     *                        \f$ \int N^T N\mathrm{d}\Omega\f$. Not used.
     * @param local_K_data    Laplacian matrix of an element, which takes the
     *         form of \f$ \int (\nabla N)^T K \nabla N\mathrm{d}\Omega\f$.
     *                        Not used.
     * @param local_b_data    Right hand side vector of an element.
     * @param local_Jac_data  Element Jacobian matrix for the Newton-Raphson
     *                        method.
     */
    void assembleWithJacobianForPressureEquations(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_xdot, const double dxdot_dx,
        const double dx_dx, std::vector<double>& local_M_data,
        std::vector<double>& local_K_data, std::vector<double>& local_b_data,
        std::vector<double>& local_Jac_data);

    unsigned getNumberOfIntegrationPoints() const override;

    typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables const&
    getMaterialStateVariablesAt(unsigned integration_point) const override;

private:
    RichardsMechanicsProcessData<DisplacementDim>& _process_data;

    std::vector<IpData, Eigen::aligned_allocator<IpData>> _ip_data;

    IntegrationMethod _integration_method;
    MeshLib::Element const& _element;
    bool const _is_axially_symmetric;
    SecondaryData<
        typename ShapeMatricesTypeDisplacement::ShapeMatrices::ShapeType>
        _secondary_data;

    static const int pressure_index = 0;
    static const int pressure_size = ShapeFunctionPressure::NPOINTS;
    static const int displacement_index = ShapeFunctionPressure::NPOINTS;
    static const int displacement_size =
        ShapeFunctionDisplacement::NPOINTS * DisplacementDim;
};

}  // namespace RichardsMechanics
}  // namespace ProcessLib

#include "RichardsMechanicsFEM-impl.h"
