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

#include <vector>

#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"

#include "ProcessLib/LIE/HydroMechanics/HydroMechanicsProcessData.h"

#include "HydroMechanicsLocalAssemblerInterface.h"
#include "IntegrationPointDataMatrix.h"

namespace ProcessLib
{
namespace LIE
{
namespace HydroMechanics
{
namespace MPL = MaterialPropertyLib;

template <typename ShapeFunctionDisplacement,
          typename ShapeFunctionPressure,
          typename IntegrationMethod,
          int GlobalDim>
class HydroMechanicsLocalAssemblerMatrix
    : public HydroMechanicsLocalAssemblerInterface
{
public:
    HydroMechanicsLocalAssemblerMatrix(
        HydroMechanicsLocalAssemblerMatrix const&) = delete;
    HydroMechanicsLocalAssemblerMatrix(HydroMechanicsLocalAssemblerMatrix&&) =
        delete;

    HydroMechanicsLocalAssemblerMatrix(
        MeshLib::Element const& e,
        std::size_t const n_variables,
        std::size_t const local_matrix_size,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        HydroMechanicsProcessData<GlobalDim>& process_data);

    void preTimestepConcrete(std::vector<double> const& /*local_x*/,
                             double const /*t*/,
                             double const /*delta_t*/) override
    {
        for (auto& data : _ip_data)
        {
            data.pushBackState();
        }
    }

protected:
    void assembleWithJacobianConcrete(double const t, double const dt,
                                      Eigen::VectorXd const& local_x,
                                      Eigen::VectorXd const& local_x_dot,
                                      Eigen::VectorXd& local_rhs,
                                      Eigen::MatrixXd& local_Jac) override;

    void assembleBlockMatricesWithJacobian(
        double const t, double const dt,
        Eigen::Ref<const Eigen::VectorXd> const& p,
        Eigen::Ref<const Eigen::VectorXd> const& p_dot,
        Eigen::Ref<const Eigen::VectorXd> const& u,
        Eigen::Ref<const Eigen::VectorXd> const& u_dot,
        Eigen::Ref<Eigen::VectorXd> rhs_p, Eigen::Ref<Eigen::VectorXd> rhs_u,
        Eigen::Ref<Eigen::MatrixXd> J_pp, Eigen::Ref<Eigen::MatrixXd> J_pu,
        Eigen::Ref<Eigen::MatrixXd> J_uu, Eigen::Ref<Eigen::MatrixXd> J_up);

    void postTimestepConcreteWithVector(
        double const t, double const dt,
        Eigen::VectorXd const& local_x) override;

    void postTimestepConcreteWithBlockVectors(
        double const t, double const dt,
        Eigen::Ref<const Eigen::VectorXd> const& p,
        Eigen::Ref<const Eigen::VectorXd> const& u);

    void setPressureOfInactiveNodes(
        double const t, Eigen::Ref<Eigen::VectorXd> p);
    void setPressureDotOfInactiveNodes(Eigen::Ref<Eigen::VectorXd> p_dot);

    // Types for displacement.
    using ShapeMatricesTypeDisplacement =
        ShapeMatrixPolicyType<ShapeFunctionDisplacement, GlobalDim>;
    using BMatricesType =
        BMatrixPolicyType<ShapeFunctionDisplacement, GlobalDim>;

    // Types for pressure.
    using ShapeMatricesTypePressure =
        ShapeMatrixPolicyType<ShapeFunctionPressure, GlobalDim>;

    using IntegrationPointDataType =
        IntegrationPointDataMatrix<BMatricesType,
                                   ShapeMatricesTypeDisplacement,
                                   ShapeMatricesTypePressure,
                                   GlobalDim,
                                   ShapeFunctionDisplacement::NPOINTS>;

    using GlobalDimVector = Eigen::Matrix<double, GlobalDim, 1>;

    HydroMechanicsProcessData<GlobalDim>& _process_data;

    std::vector<IntegrationPointDataType,
                Eigen::aligned_allocator<IntegrationPointDataType>>
        _ip_data;

    static const int pressure_index = 0;
    static const int pressure_size = ShapeFunctionPressure::NPOINTS;
    static const int displacement_index = ShapeFunctionPressure::NPOINTS;
    static const int displacement_size =
        ShapeFunctionDisplacement::NPOINTS * GlobalDim;
    static const int kelvin_vector_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(GlobalDim);
};

}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib

#include "HydroMechanicsLocalAssemblerMatrix-impl.h"
