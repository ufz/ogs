/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "HydroMechanicsLocalAssemblerInterface.h"
#include "IntegrationPointDataFracture.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/LIE/Common/HMatrixUtils.h"
#include "ProcessLib/LIE/HydroMechanics/HydroMechanicsProcessData.h"

namespace ProcessLib
{
namespace LIE
{
namespace HydroMechanics
{
template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,

          int DisplacementDim>
class HydroMechanicsLocalAssemblerFracture
    : public HydroMechanicsLocalAssemblerInterface
{
public:
    HydroMechanicsLocalAssemblerFracture(
        HydroMechanicsLocalAssemblerFracture const&) = delete;
    HydroMechanicsLocalAssemblerFracture(
        HydroMechanicsLocalAssemblerFracture&&) = delete;

    HydroMechanicsLocalAssemblerFracture(
        MeshLib::Element const& e,
        std::size_t const local_matrix_size,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        HydroMechanicsProcessData<DisplacementDim>& process_data);

    void preTimestepConcrete(std::vector<double> const& /*local_x*/,
                             double const /*t*/,
                             double const /*delta_t*/) override
    {
        for (auto& data : _ip_data)
        {
            data.pushBackState();
        }
    }

    void postTimestepConcreteWithVector(
        const double t, double const dt,
        Eigen::VectorXd const& local_x) override;

    std::vector<double> const& getIntPtSigma(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        cache.resize(0);
        return cache;
    }

    std::vector<double> const& getIntPtEpsilon(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        cache.resize(0);
        return cache;
    }

    std::vector<double> const& getIntPtDarcyVelocity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        cache.resize(0);
        return cache;
    }

    std::vector<double> const& getIntPtFractureVelocity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override;

    std::vector<double> const& getIntPtFractureStress(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override;

    std::vector<double> const& getIntPtFractureAperture(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override;

    std::vector<double> const& getIntPtFracturePermeability(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override;

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _secondary_data.N[integration_point];

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

private:
    void assembleWithJacobianConcrete(double const t, double const dt,
                                      Eigen::VectorXd const& local_x,
                                      Eigen::VectorXd const& local_x_prev,
                                      Eigen::VectorXd& local_b,
                                      Eigen::MatrixXd& local_J) override;

    void assembleBlockMatricesWithJacobian(
        double const t, double const dt,
        Eigen::Ref<const Eigen::VectorXd> const& p,
        Eigen::Ref<const Eigen::VectorXd> const& p_prev,
        Eigen::Ref<const Eigen::VectorXd> const& g,
        Eigen::Ref<const Eigen::VectorXd> const& g_prev,
        Eigen::Ref<Eigen::VectorXd> rhs_p, Eigen::Ref<Eigen::VectorXd> rhs_g,
        Eigen::Ref<Eigen::MatrixXd> J_pp, Eigen::Ref<Eigen::MatrixXd> J_pg,
        Eigen::Ref<Eigen::MatrixXd> J_gg, Eigen::Ref<Eigen::MatrixXd> J_gp);

    // Types for displacement.
    using ShapeMatricesTypeDisplacement =
        ShapeMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;
    using HMatricesType =
        HMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;
    using HMatrixType = typename HMatricesType::HMatrixType;

    // Types for pressure.
    using ShapeMatricesTypePressure =
        ShapeMatrixPolicyType<ShapeFunctionPressure, DisplacementDim>;

    // Types for the integration point data
    using IntegrationPointDataType = IntegrationPointDataFracture<
        HMatricesType, ShapeMatricesTypeDisplacement, ShapeMatricesTypePressure,
        DisplacementDim>;

    using GlobalDimVectorType = Eigen::Matrix<double, DisplacementDim, 1>;

private:
    HydroMechanicsProcessData<DisplacementDim>& _process_data;
    FractureProperty const* _fracture_property = nullptr;

    std::vector<IntegrationPointDataType,
                Eigen::aligned_allocator<IntegrationPointDataType>>
        _ip_data;

    static const int pressure_index = 0;
    static const int pressure_size = ShapeFunctionPressure::NPOINTS;
    static const int displacement_index = ShapeFunctionPressure::NPOINTS;
    static const int displacement_size =
        ShapeFunctionDisplacement::NPOINTS * DisplacementDim;

    SecondaryData<
        typename ShapeMatricesTypeDisplacement::ShapeMatrices::ShapeType>
        _secondary_data;
};

}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib

#include "HydroMechanicsLocalAssemblerFracture-impl.h"
