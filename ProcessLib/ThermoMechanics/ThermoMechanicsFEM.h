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

#include <memory>
#include <vector>

#include "LocalAssemblerInterface.h"
#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "MathLib/KelvinVector.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/Integration/GenericIntegrationMethod.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ThermoMechanicsProcessData.h"

namespace ProcessLib
{
namespace ThermoMechanics
{
namespace MPL = MaterialPropertyLib;

template <typename BMatricesType, typename ShapeMatricesType,
          int DisplacementDim>
struct IntegrationPointData final
{
    explicit IntegrationPointData(
        MaterialLib::Solids::MechanicsBase<DisplacementDim> const&
            solid_material)
        : solid_material(solid_material),
          material_state_variables(
              solid_material.createMaterialStateVariables())
    {
    }

    // total stress
    typename BMatricesType::KelvinVectorType sigma, sigma_prev;

    // total strain
    typename BMatricesType::KelvinVectorType eps, eps_prev;

    // mechanical strain
    typename BMatricesType::KelvinVectorType eps_m, eps_m_prev;

    MaterialLib::Solids::MechanicsBase<DisplacementDim> const& solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        material_state_variables;

    double integration_weight;
    typename ShapeMatricesType::NodalRowVectorType N;
    typename ShapeMatricesType::GlobalDimNodalMatrixType dNdx;

    void pushBackState()
    {
        eps_prev = eps;
        eps_m_prev = eps_m;
        sigma_prev = sigma;
        material_state_variables->pushBackState();
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

/// Used for the extrapolation of the integration point values. It is ordered
/// (and stored) by integration points.
template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType, Eigen::aligned_allocator<ShapeMatrixType>> N;
};

/**
 * \brief Local assembler of ThermoMechanics process.
 *
 *  In this assembler, the solid density can be expressed as an exponential
 *  function of temperature, if a temperature dependent solid density is
 *  assumed.
 *  The theory behind this is given below.
 *
 *  During the temperature variation, the solid mass balance is
 *  \f[
 *    \rho_s^{n+1}  V^{n+1} = \rho_s^n  V^n,
 *  \f]
 *  with \f$\rho_s\f$ the density,  \f$V\f$, the volume, and \f$n\f$ the index
 *  of the time step.
 *  Under pure thermo-mechanics condition, the volume change along with
 *  the temperature change is given by
 *  \f[
 *   V^{n+1} = V^{n} +  \alpha_T  dT  V^{n} = (1+\alpha_T  dT)V^{n},
 *  \f]
 *  where \f$\alpha_T\f$ is the volumetric thermal expansivity of the solid.
 *  This gives
 *  \f[ \rho_s^{n+1} +  \alpha_T  dT \rho_s^{n+1} =  \rho_s^n. \f]
 *   Therefore, we obtain the differential expression of the
 *   temperature dependent solid density
 *   as
 *   \f[
 *        \frac{d \rho_s}{d T} = -\alpha_T \rho_s.
 *   \f]
 *    The solution of the above ODE, i.e the density expression, is given by
 *   \f[
 *        \rho_s = {\rho_0} \mathrm{exp} (- \alpha_T (T-T0)),
 *   \f]
 *   with reference density \f$\rho_0\f$ at a reference temperature of
 *  \f$T_0\f$.
 *
 *   An MPL property with the type of **Exponential** (see
 *   MaterialPropertyLib::Exponential) can be used for the
 *   temperature dependent solid property.
 */

template <typename ShapeFunction, int DisplacementDim>
class ThermoMechanicsLocalAssembler
    : public ThermoMechanicsLocalAssemblerInterface<DisplacementDim>
{
public:
    using ShapeMatricesType =
        ShapeMatrixPolicyType<ShapeFunction, DisplacementDim>;

    // Types for displacement.
    // (Higher order elements = ShapeFunction).
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;
    using GlobalDimMatrixType = typename ShapeMatricesType::GlobalDimMatrixType;
    using BMatricesType = BMatrixPolicyType<ShapeFunction, DisplacementDim>;

    static int const KelvinVectorSize =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;

    using NodalForceVectorType = typename BMatricesType::NodalForceVectorType;
    using RhsVector = typename ShapeMatricesType::template VectorType<
        ShapeFunction::NPOINTS + ShapeFunction::NPOINTS * DisplacementDim>;
    using JacobianMatrix = typename ShapeMatricesType::template MatrixType<
        ShapeFunction::NPOINTS + ShapeFunction::NPOINTS * DisplacementDim,
        ShapeFunction::NPOINTS + ShapeFunction::NPOINTS * DisplacementDim>;
    using IpData =
        IntegrationPointData<BMatricesType, ShapeMatricesType, DisplacementDim>;

    ThermoMechanicsLocalAssembler(ThermoMechanicsLocalAssembler const&) =
        delete;
    ThermoMechanicsLocalAssembler(ThermoMechanicsLocalAssembler&&) = delete;

    ThermoMechanicsLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        ThermoMechanicsProcessData<DisplacementDim>& process_data);

    /// Returns number of read integration points.
    std::size_t setIPDataInitialConditions(
        std::string const& name,
        double const* values,
        int const integration_order) override;

    void assemble(double const /*t*/, double const /*dt*/,
                  std::vector<double> const& /*local_x*/,
                  std::vector<double> const& /*local_x_prev*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_rhs_data*/) override
    {
        OGS_FATAL(
            "ThermoMechanicsLocalAssembler: assembly without jacobian is not "
            "implemented.");
    }

    void assembleWithJacobianForStaggeredScheme(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_x_prev, int const process_id,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data,
        std::vector<double>& local_Jac_data) override;

    void assembleWithJacobian(double const t, double const dt,
                              std::vector<double> const& local_x,
                              std::vector<double> const& local_x_prev,
                              std::vector<double>& local_M_data,
                              std::vector<double>& local_K_data,
                              std::vector<double>& local_rhs_data,
                              std::vector<double>& local_Jac_data) override;

    void initializeConcrete() override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto& ip_data = _ip_data[ip];

            ParameterLib::SpatialPosition const x_position{
                std::nullopt, _element.getID(), ip,
                MathLib::Point3d(
                    NumLib::interpolateCoordinates<ShapeFunction,
                                                   ShapeMatricesType>(
                        _element, ip_data.N))};

            /// Set initial stress from parameter.
            if (_process_data.initial_stress != nullptr)
            {
                ip_data.sigma =
                    MathLib::KelvinVector::symmetricTensorToKelvinVector<
                        DisplacementDim>((*_process_data.initial_stress)(
                        std::numeric_limits<
                            double>::quiet_NaN() /* time independent */,
                        x_position));
            }

            double const t = 0;  // TODO (naumov) pass t from top
            ip_data.solid_material.initializeInternalStateVariables(
                t, x_position, *ip_data.material_state_variables);

            ip_data.pushBackState();
        }
    }

    void postTimestepConcrete(Eigen::VectorXd const& /*local_x*/,
                              Eigen::VectorXd const& /*local_x_prev*/,
                              double const /*t*/, double const /*dt*/,
                              bool const /*use_monolithic_scheme*/,
                              int const /*process_id*/) override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data[ip].pushBackState();
        }
    }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _secondary_data.N[integration_point];

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

private:
    /**
     * Assemble local matrices and vectors arise from the linearized discretized
     * weak form of the residual of the momentum balance equation,
     *      \f[
     *            \nabla (\sigma - \mathbf{D} \alpha_T (T-T_0) \mathrm{I}) = f
     *      \f]
     * where \f$ \sigma\f$ is the effective stress tensor, \f$\mathbf{D}\f$ is
     * the tangential operator, \f$T\f$ is the  temperature, \f$T_0\f$ is the
     * initial temperature, \f$\alpha_T\f$ is the linear thermal expansion,
     * \f$\mathrm{I}\f$ is the identity tensor, and \f$f\f$ is the body force.
     *
     * @param t               Time
     * @param dt              Time increment
     * @param local_x         Nodal values of \f$x\f$ of all processes of an
     * element.
     * @param local_x_prev    Nodal values of \f$x_{prev}\f$ of all processes of
     * an element.
     * @param local_b_data    Right hand side vector of an element.
     * @param local_Jac_data  Element Jacobian matrix for the Newton-Raphson
     *                        method.
     */

    void assembleWithJacobianForDeformationEquations(
        const double t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_x_prev, std::vector<double>& local_b_data,
        std::vector<double>& local_Jac_data);

    /**
     * Assemble local matrices and vectors arise from the linearized discretized
     * weak form of the residual of the energy balance equation,
     *      \f[
     *          \rho c_p \cdot{T} - \nabla (\mathbf{K} (\nabla T) = Q_T
     *      \f]
     * where \f$ rho\f$ is the solid density, \f$ c_p\f$ is the specific heat
     * capacity, \f$ \mathbf{K} \f$ is the thermal conductivity, and \f$ Q_T\f$
     * is the source/sink term.
     *
     * @param t               Time
     * @param dt              Time increment
     * @param local_x         Nodal values of \f$x\f$ of all processes of an
     * element.
     * @param local_x_prev    Nodal values of \f$x_{prev}\f$ of all processes of
     * an element.
     * @param local_b_data    Right hand side vector of an element.
     * @param local_Jac_data  Element Jacobian matrix for the Newton-Raphson
     *                        method.
     */
    void assembleWithJacobianForHeatConductionEquations(
        const double t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_x_prev, std::vector<double>& local_b_data,
        std::vector<double>& local_Jac_data);

    std::size_t setSigma(double const* values);

    // TODO (naumov) This method is same as getIntPtSigma but for arguments and
    // the ordering of the cache_mat.
    // There should be only one.
    std::vector<double> getSigma() const override;

    std::vector<double> const& getIntPtSigma(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override;

    std::size_t setEpsilon(double const* values);

    std::vector<double> getEpsilon() const override;

    std::vector<double> const& getIntPtEpsilon(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override;

    std::vector<double> const& getIntPtEpsilonMechanical(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override;

    std::size_t setEpsilonMechanical(double const* values);

    std::vector<double> getEpsilonMechanical() const override;

    unsigned getNumberOfIntegrationPoints() const override;

    int getMaterialID() const override;

    typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables const&
    getMaterialStateVariablesAt(unsigned integration_point) const override;

private:
    ThermoMechanicsProcessData<DisplacementDim>& _process_data;

    std::vector<IpData, Eigen::aligned_allocator<IpData>> _ip_data;

    NumLib::GenericIntegrationMethod const& _integration_method;
    MeshLib::Element const& _element;
    SecondaryData<typename ShapeMatrices::ShapeType> _secondary_data;
    bool const _is_axially_symmetric;

    static const int temperature_index = 0;
    static const int temperature_size = ShapeFunction::NPOINTS;
    static const int displacement_index = ShapeFunction::NPOINTS;
    static const int displacement_size =
        ShapeFunction::NPOINTS * DisplacementDim;
};

}  // namespace ThermoMechanics
}  // namespace ProcessLib

#include "ThermoMechanicsFEM-impl.h"
