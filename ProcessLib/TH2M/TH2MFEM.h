/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <vector>

#include "ConstitutiveVariables.h"
#include "IntegrationPointData.h"
#include "LocalAssemblerInterface.h"
#include "MaterialLib/PhysicalConstant.h"
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
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"
#include "ProcessLib/Utils/TransposeInPlace.h"
#include "TH2MProcessData.h"

namespace ProcessLib
{
namespace TH2M
{
/// Used by for extrapolation of the integration point values. It is ordered
/// (and stored) by integration points.
template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType, Eigen::aligned_allocator<ShapeMatrixType>> N_u;
};

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
class TH2MLocalAssembler : public LocalAssemblerInterface
{
public:
    using ShapeMatricesTypeDisplacement =
        ShapeMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;

    using ShapeMatricesTypePressure =
        ShapeMatrixPolicyType<ShapeFunctionPressure, DisplacementDim>;

    template <int N>
    using VectorType =
        typename ShapeMatricesTypePressure::template VectorType<N>;

    template <int M, int N>
    using MatrixType =
        typename ShapeMatricesTypePressure::template MatrixType<M, N>;

    using GlobalDimMatrixType =
        typename ShapeMatricesTypePressure::GlobalDimMatrixType;

    using GlobalDimVectorType =
        typename ShapeMatricesTypePressure::GlobalDimVectorType;

    static int const KelvinVectorSize =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using SymmetricTensor = Eigen::Matrix<double, KelvinVectorSize, 1>;

    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;

    TH2MLocalAssembler(TH2MLocalAssembler const&) = delete;
    TH2MLocalAssembler(TH2MLocalAssembler&&) = delete;

    TH2MLocalAssembler(MeshLib::Element const& e,
                       std::size_t const /*local_matrix_size*/,
                       bool const is_axially_symmetric,
                       unsigned const integration_order,
                       TH2MProcessData<DisplacementDim>& process_data);

private:
    /// \return the number of read integration points.
    std::size_t setIPDataInitialConditions(
        std::string const& name,
        double const* values,
        int const integration_order) override;

    void setInitialConditionsConcrete(std::vector<double> const& local_x,
                                      double const t,
                                      bool const /*use_monolithic_scheme*/,
                                      int const /*process_id*/) override;

    void assemble(double const /*t*/, double const /*dt*/,
                  std::vector<double> const& /*local_x*/,
                  std::vector<double> const& /*local_xdot*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_rhs_data*/) override;

    void assembleWithJacobian(double const t, double const dt,
                              std::vector<double> const& local_x,
                              std::vector<double> const& local_xdot,
                              std::vector<double>& /*local_M_data*/,
                              std::vector<double>& /*local_K_data*/,
                              std::vector<double>& local_rhs_data,
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

    std::vector<double> const& getIntPtDarcyVelocityGas(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override;
    std::vector<double> const& getIntPtDarcyVelocityLiquid(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override;

    std::vector<ConstitutiveVariables<DisplacementDim>>
    updateConstitutiveVariables(Eigen::VectorXd const& local_x,
                                Eigen::VectorXd const& local_x_dot,
                                double const t, double const dt);

    std::size_t setSigma(double const* values)
    {
        return ProcessLib::setIntegrationPointKelvinVectorData<DisplacementDim>(
            values, _ip_data, &IpData::sigma_eff);
    }

    // TODO (naumov) This method is same as getIntPtSigma but for arguments and
    // the ordering of the cache_mat.
    // There should be only one.
    std::vector<double> getSigma() const override
    {
        {
            constexpr int kelvin_vector_size =
                MathLib::KelvinVector::kelvin_vector_dimensions(
                    DisplacementDim);

            return transposeInPlace<kelvin_vector_size>(
                [this](std::vector<double>& values)
                { return getIntPtSigma(0, {}, {}, values); });
        }
    }

    std::vector<double> const& getIntPtSigma(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
            _ip_data, &IpData::sigma_eff, cache);
    }

    std::vector<double> getEpsilon() const override
    {
        constexpr int kelvin_vector_size =
            MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);

        return transposeInPlace<kelvin_vector_size>(
            [this](std::vector<double>& values)
            { return getIntPtEpsilon(0, {}, {}, values); });
    }

    virtual std::vector<double> const& getIntPtEpsilon(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
            _ip_data, &IpData::eps, cache);
    }

    // TODO: Here is some refactoring potential. All secondary variables could
    // be stored in some container to avoid defining one method for each
    // variable.

    virtual std::vector<double> const& getIntPtLiquidDensity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(_ip_data,
                                                         &IpData::rhoLR, cache);
    }

    virtual std::vector<double> const& getIntPtGasDensity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(_ip_data,
                                                         &IpData::rhoGR, cache);
    }

    virtual std::vector<double> const& getIntPtSolidDensity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(_ip_data,
                                                         &IpData::rhoSR, cache);
    }

    virtual std::vector<double> const& getIntPtVapourPressure(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(_ip_data,
                                                         &IpData::pWGR, cache);
    }

    virtual std::vector<double> const& getIntPtPorosity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(_ip_data, &IpData::phi,
                                                         cache);
    }

    std::vector<double> getSaturation() const override
    {
        std::vector<double> result;
        getIntPtSaturation(0, {}, {}, result);
        return result;
    }

    virtual std::vector<double> const& getIntPtSaturation(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(_ip_data, &IpData::s_L,
                                                         cache);
    }

    virtual std::vector<double> const& getIntPtMoleFractionGas(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(_ip_data,
                                                         &IpData::xnCG, cache);
    }
    virtual std::vector<double> const& getIntPtMassFractionGas(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(_ip_data,
                                                         &IpData::xmCG, cache);
    }
    virtual std::vector<double> const& getIntPtMassFractionLiquid(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(_ip_data,
                                                         &IpData::xmWL, cache);
    }

    virtual std::vector<double> const& getIntPtRelativePermeabilityGas(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(
            _ip_data, &IpData::k_rel_G, cache);
    }
    virtual std::vector<double> const& getIntPtRelativePermeabilityLiquid(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(
            _ip_data, &IpData::k_rel_L, cache);
    }

    virtual std::vector<double> const& getIntPtEnthalpyGas(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(_ip_data, &IpData::h_G,
                                                         cache);
    }
    virtual std::vector<double> const& getIntPtEnthalpyLiquid(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(_ip_data, &IpData::h_L,
                                                         cache);
    }
    virtual std::vector<double> const& getIntPtEnthalpySolid(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(_ip_data, &IpData::h_S,
                                                         cache);
    }

private:
    TH2MProcessData<DisplacementDim>& _process_data;

    using BMatricesType =
        BMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;
    using IpData =
        IntegrationPointData<BMatricesType, ShapeMatricesTypeDisplacement,
                             ShapeMatricesTypePressure, DisplacementDim,
                             ShapeFunctionDisplacement::NPOINTS>;
    std::vector<IpData, Eigen::aligned_allocator<IpData>> _ip_data;

    IntegrationMethod _integration_method;
    MeshLib::Element const& _element;
    bool const _is_axially_symmetric;
    SecondaryData<
        typename ShapeMatricesTypeDisplacement::ShapeMatrices::ShapeType>
        _secondary_data;

    // The shape function of pressure has the same form with the shape function
    // of temperature
    static const int gas_pressure_index = 0;
    static const int gas_pressure_size = ShapeFunctionPressure::NPOINTS;
    static const int capillary_pressure_index = ShapeFunctionPressure::NPOINTS;
    static const int capillary_pressure_size = ShapeFunctionPressure::NPOINTS;
    static const int temperature_index = 2 * ShapeFunctionPressure::NPOINTS;
    static const int temperature_size = ShapeFunctionPressure::NPOINTS;
    static const int displacement_index = ShapeFunctionPressure::NPOINTS * 3;
    static const int displacement_size =
        ShapeFunctionDisplacement::NPOINTS * DisplacementDim;

    static const int C_index = 0;
    static const int C_size = ShapeFunctionPressure::NPOINTS;
    static const int W_index = ShapeFunctionPressure::NPOINTS;
    static const int W_size = ShapeFunctionPressure::NPOINTS;
};

}  // namespace TH2M
}  // namespace ProcessLib

#include "TH2MFEM-impl.h"
