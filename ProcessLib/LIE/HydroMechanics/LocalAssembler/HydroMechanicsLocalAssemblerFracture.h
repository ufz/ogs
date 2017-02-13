/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "NumLib/Fem/ShapeMatrixPolicy.h"

#include "ProcessLib/LIE/Common/HMatrixUtils.h"
#include "ProcessLib/LIE/HydroMechanics/HydroMechanicsProcessData.h"

#include "IntegrationPointDataFracture.h"
#include "HydroMechanicsLocalAssemblerInterface.h"

namespace ProcessLib
{
namespace LIE
{
namespace HydroMechanics
{

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, unsigned GlobalDim>
class HydroMechanicsLocalAssemblerFracture
    : public HydroMechanicsLocalAssemblerInterface
{
public:
    HydroMechanicsLocalAssemblerFracture(HydroMechanicsLocalAssemblerFracture const&) =
        delete;
    HydroMechanicsLocalAssemblerFracture(HydroMechanicsLocalAssemblerFracture&&) = delete;

    HydroMechanicsLocalAssemblerFracture(
        MeshLib::Element const& e,
        std::size_t const local_matrix_size,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        HydroMechanicsProcessData<GlobalDim>& process_data);

    void preTimestepConcrete(std::vector<double> const& /*local_x*/,
                             double const /*t*/,
                             double const /*delta_t*/) override
    {
        for (auto &data : _ip_data)
            data.pushBackState();
    }

    void computeSecondaryVariableConcreteWithVector(
        const double t, Eigen::VectorXd const& local_x) override;

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _ip_data[integration_point].N_p;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

private:
    void assembleWithJacobianConcrete(
        double const t,
        Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_xdot,
        Eigen::VectorXd& local_b,
        Eigen::MatrixXd& local_J) override;

    void assembleBlockMatricesWithJacobian(
        double const t,
        Eigen::Ref<const Eigen::VectorXd> const& p,
        Eigen::Ref<const Eigen::VectorXd> const& p_dot,
        Eigen::Ref<const Eigen::VectorXd> const& g,
        Eigen::Ref<const Eigen::VectorXd> const& g_dot,
        Eigen::Ref<Eigen::VectorXd> rhs_p,
        Eigen::Ref<Eigen::VectorXd> rhs_g,
        Eigen::Ref<Eigen::MatrixXd> J_pp,
        Eigen::Ref<Eigen::MatrixXd> J_pg,
        Eigen::Ref<Eigen::MatrixXd> J_gg,
        Eigen::Ref<Eigen::MatrixXd> J_gp);

    // Types for displacement.
    using ShapeMatricesTypeDisplacement =
        ShapeMatrixPolicyType<ShapeFunctionDisplacement, GlobalDim>;
    using HMatricesType = HMatrixPolicyType<ShapeFunctionDisplacement, GlobalDim>;
    using HMatrixType = typename HMatricesType::HMatrixType;

    // Types for pressure.
    using ShapeMatricesTypePressure =
        ShapeMatrixPolicyType<ShapeFunctionPressure, GlobalDim>;

    // Types for the integration point data
    using IntegrationPointDataType =
        IntegrationPointDataFracture<HMatricesType,
                                     ShapeMatricesTypeDisplacement,
                                     ShapeMatricesTypePressure, GlobalDim>;

    HydroMechanicsProcessData<GlobalDim>& _process_data;

    std::vector<IntegrationPointDataType,
        Eigen::aligned_allocator<IntegrationPointDataType>> _ip_data;

    static const int pressure_index = 0;
    static const int pressure_size = ShapeFunctionPressure::NPOINTS;
    static const int displacement_index = ShapeFunctionPressure::NPOINTS;
    static const int displacement_size =
        ShapeFunctionDisplacement::NPOINTS * GlobalDim;
};

}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib

#include "HydroMechanicsLocalAssemblerFracture-impl.h"
