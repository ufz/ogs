/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <vector>

#include "ConstitutiveRelations/Base.h"
#include "ConstitutiveRelations/ConstitutiveData.h"
#include "IntegrationPointData.h"
#include "LocalAssemblerInterface.h"
#include "MaterialLib/MPL/VariableType.h"
#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MathLib/EigenBlockMatrixView.h"
#include "MathLib/KelvinVector.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/Integration/GenericIntegrationMethod.h"
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
          int DisplacementDim>
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

    static constexpr auto& N_u_op = MathLib::eigenBlockMatrixView<
        DisplacementDim,
        typename ShapeMatricesTypeDisplacement::NodalRowVectorType>;

    RichardsMechanicsLocalAssembler(RichardsMechanicsLocalAssembler const&) =
        delete;
    RichardsMechanicsLocalAssembler(RichardsMechanicsLocalAssembler&&) = delete;

    RichardsMechanicsLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        RichardsMechanicsProcessData<DisplacementDim>& process_data);

    /// \return the number of read integration points.
    std::size_t setIPDataInitialConditions(
        std::string_view const name,
        double const* values,
        int const integration_order) override;

    void setInitialConditionsConcrete(Eigen::VectorXd const local_x,
                                      double const t,
                                      int const process_id) override;

    void assemble(double const t, double const dt,
                  std::vector<double> const& local_x,
                  std::vector<double> const& local_x_prev,
                  std::vector<double>& local_M_data,
                  std::vector<double>& local_K_data,
                  std::vector<double>& local_rhs_data) override;

    void assembleWithJacobian(double const t, double const dt,
                              std::vector<double> const& local_x,
                              std::vector<double> const& local_x_prev,
                              std::vector<double>& /*local_M_data*/,
                              std::vector<double>& /*local_K_data*/,
                              std::vector<double>& local_rhs_data,
                              std::vector<double>& local_Jac_data) override;

    void assembleWithJacobianForStaggeredScheme(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_x_prev, int const process_id,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data,
        std::vector<double>& local_Jac_data) override;

    void initializeConcrete() override
    {
        unsigned const n_integration_points =
            this->integration_method_.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto& SD = this->current_states_[ip];
            auto& ip_data = _ip_data[ip];

            ParameterLib::SpatialPosition const x_position{
                std::nullopt, this->element_.getID(), ip,
                MathLib::Point3d(NumLib::interpolateCoordinates<
                                 ShapeFunctionDisplacement,
                                 ShapeMatricesTypeDisplacement>(this->element_,
                                                                ip_data.N_u))};

            /// Set initial stress from parameter.
            if (this->process_data_.initial_stress != nullptr)
            {
                std::get<ProcessLib::ThermoRichardsMechanics::
                             ConstitutiveStress_StrainTemperature::
                                 EffectiveStressData<DisplacementDim>>(SD)
                    .sigma_eff =
                    MathLib::KelvinVector::symmetricTensorToKelvinVector<
                        DisplacementDim>((*this->process_data_.initial_stress)(
                        std::numeric_limits<
                            double>::quiet_NaN() /* time independent */,
                        x_position));
            }

            double const t = 0;  // TODO (naumov) pass t from top
            this->solid_material_.initializeInternalStateVariables(
                t, x_position,
                *this->material_states_[ip].material_state_variables);

            this->material_states_[ip].pushBackState();

            this->prev_states_[ip] = SD;
        }
    }

    void postTimestepConcrete(Eigen::VectorXd const& /*local_x*/,
                              Eigen::VectorXd const& /*local_x_prev*/,
                              double const /*t*/, double const /*dt*/,
                              int const /*process_id*/) override
    {
        unsigned const n_integration_points =
            this->integration_method_.getNumberOfPoints();

        for (auto& s : this->material_states_)
        {
            s.pushBackState();
        }

        // TODO move to the local assembler interface
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            this->prev_states_[ip] = this->current_states_[ip];
        }
    }

    void computeSecondaryVariableConcrete(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_x_prev) override;

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N_u = _secondary_data.N_u[integration_point];

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N_u.data(), N_u.size());
    }

    int getMaterialID() const override;

    std::vector<double> getMaterialStateVariableInternalState(
        std::function<std::span<double>(
            typename MaterialLib::Solids::MechanicsBase<DisplacementDim>::
                MaterialStateVariables&)> const& get_values_span,
        int const& n_components) const override;

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
     * @param local_x_prev    Nodal values of \f$x_{prev}\f$ of an element.
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
        Eigen::VectorXd const& local_x_prev, std::vector<double>& local_M_data,
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
     * @param local_x_prev    Nodal values of \f$x_{prev}\f$ of an element.
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
        Eigen::VectorXd const& local_x_prev, std::vector<double>& local_M_data,
        std::vector<double>& local_K_data, std::vector<double>& local_b_data,
        std::vector<double>& local_Jac_data);

    unsigned getNumberOfIntegrationPoints() const override;

    typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables const&
    getMaterialStateVariablesAt(unsigned integration_point) const override;

private:
    static void assembleWithJacobianEvalConstitutiveSetting(
        double const t, double const dt,
        ParameterLib::SpatialPosition const& x_position, IpData& ip_data,
        MPL::VariableArray& variables, MPL::VariableArray& variables_prev,
        MPL::Medium const* const medium, TemperatureData const T_data,
        CapillaryPressureData<DisplacementDim> const& p_cap_data,
        ConstitutiveData<DisplacementDim>& CD,
        StatefulData<DisplacementDim>& SD,
        StatefulDataPrev<DisplacementDim> const& SD_prev,
        std::optional<MicroPorosityParameters> const& micro_porosity_parameters,
        MaterialLib::Solids::MechanicsBase<DisplacementDim> const&
            solid_material,
        ProcessLib::ThermoRichardsMechanics::MaterialStateData<DisplacementDim>&
            material_state_data);

    std::vector<IpData, Eigen::aligned_allocator<IpData>> _ip_data;

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
