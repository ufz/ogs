/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file
 *
 * Created on August 19, 2016, 2:28 PM
 */

#pragma once

#include <map>
#include <unordered_map>
#include <vector>
#include <typeindex>

#include "MaterialLib/PorousMedium/Permeability/Permeability.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "LiquidFlowMaterialProperties.h"

namespace ProcessLib
{
namespace LiquidFlow
{
template <typename NodalRowVectorType, typename GlobalDimNodalMatrixType>
struct IntegrationPointData final
{
    explicit IntegrationPointData(NodalRowVectorType const& N_,
                                  GlobalDimNodalMatrixType const& dNdx_,
                                  double const& integration_weight_)
        : N(N_),
          dNdx(dNdx_),
          integration_weight(integration_weight_)
    {}

    NodalRowVectorType const N;
    GlobalDimNodalMatrixType const dNdx;
    double const integration_weight;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

const unsigned NUM_NODAL_DOF = 1;

class LiquidFlowLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    virtual std::vector<double> const& getIntPtDarcyVelocity(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const = 0;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class LiquidFlowLocalAssembler : public LiquidFlowLocalAssemblerInterface
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalAssemblerTraits = ProcessLib::LocalAssemblerTraits<
        ShapeMatricesType, ShapeFunction::NPOINTS, NUM_NODAL_DOF, GlobalDim>;

    using NodalMatrixType = typename LocalAssemblerTraits::LocalMatrix;
    using NodalVectorType = typename LocalAssemblerTraits::LocalVector;
    using NodalRowVectorType = typename ShapeMatricesType::NodalRowVectorType;
    using GlobalDimNodalMatrixType =
        typename ShapeMatricesType::GlobalDimNodalMatrixType;

    using MatrixOfVelocityAtIntegrationPoints = Eigen::Map<
        Eigen::Matrix<double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>;

public:
    LiquidFlowLocalAssembler(
        MeshLib::Element const& element,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        int const gravitational_axis_id,
        double const gravitational_acceleration,
        double const reference_temperature,
        LiquidFlowMaterialProperties const& material_propertries)
        : _element(element),
          _integration_method(integration_order),
          _gravitational_axis_id(gravitational_axis_id),
          _gravitational_acceleration(gravitational_acceleration),
          _reference_temperature(reference_temperature),
          _material_properties(material_propertries)
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        _ip_data.reserve(n_integration_points);

        auto const& shape_matrices =
            initShapeMatrices<ShapeFunction, ShapeMatricesType,
                              IntegrationMethod, GlobalDim>(
                element, is_axially_symmetric, _integration_method);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data.emplace_back(
                shape_matrices[ip].N, shape_matrices[ip].dNdx,
                _integration_method.getWeightedPoint(ip).getWeight() *
                    shape_matrices[ip].integralMeasure *
                    shape_matrices[ip].detJ);
        }
    }

    void assemble(double const t, std::vector<double> const& local_x,
                  std::vector<double>& local_M_data,
                  std::vector<double>& local_K_data,
                  std::vector<double>& local_b_data) override;


    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _ip_data[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    std::vector<double> const& getIntPtDarcyVelocity(
        const double t,
        GlobalVector const& current_solution,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        std::vector<double>& velocity_cache) const override;

private:
    MeshLib::Element const& _element;

    IntegrationMethod const _integration_method;
    std::vector<
        IntegrationPointData<NodalRowVectorType, GlobalDimNodalMatrixType>,
        Eigen::aligned_allocator<
            IntegrationPointData<NodalRowVectorType, GlobalDimNodalMatrixType>>>
        _ip_data;

    /**
     *  Calculator of the Laplacian and the gravity term for anisotropic
     *  permeability tensor
     */
    struct AnisotropicCalculator
    {
        static void calculateLaplacianAndGravityTerm(
            Eigen::Map<NodalMatrixType>& local_K,
            Eigen::Map<NodalVectorType>& local_b,
            IntegrationPointData<NodalRowVectorType,
                                 GlobalDimNodalMatrixType> const& ip_data,
            Eigen::MatrixXd const& permeability, double const mu,
            double const rho_g, int const gravitational_axis_id);

        static void calculateVelocity(
            unsigned const ip, Eigen::Map<const NodalVectorType> const& local_p,
            IntegrationPointData<NodalRowVectorType,
                                 GlobalDimNodalMatrixType> const& ip_data,
            Eigen::MatrixXd const& permeability, double const mu,
            double const rho_g, int const gravitational_axis_id,
            MatrixOfVelocityAtIntegrationPoints& darcy_velocity_at_ips);
    };

    /**
     *  Calculator of the Laplacian and the gravity term for isotropic
     *  permeability tensor
     */
    struct IsotropicCalculator
    {
        static void calculateLaplacianAndGravityTerm(
            Eigen::Map<NodalMatrixType>& local_K,
            Eigen::Map<NodalVectorType>& local_b,
            IntegrationPointData<NodalRowVectorType,
                                 GlobalDimNodalMatrixType> const& ip_data,
            Eigen::MatrixXd const& permeability, double const mu,
            double const rho_g, int const gravitational_axis_id);

        static void calculateVelocity(
            unsigned const ip, Eigen::Map<const NodalVectorType> const& local_p,
            IntegrationPointData<NodalRowVectorType,
                                 GlobalDimNodalMatrixType> const& ip_data,
            Eigen::MatrixXd const& permeability, double const mu,
            double const rho_g, int const gravitational_axis_id,
            MatrixOfVelocityAtIntegrationPoints& darcy_velocity_at_ips);
    };

    template <typename LaplacianGravityVelocityCalculator>
    void assembleMatrixAndVector(const int material_id, double const t,
                                 std::vector<double> const& local_x,
                                 std::vector<double>& local_M_data,
                                 std::vector<double>& local_K_data,
                                 std::vector<double>& local_b_data);

    template <typename LaplacianGravityVelocityCalculator>
    void computeDarcyVelocityLocal(
        const int material_id,
        const double t,
        std::vector<double> const& local_x,
        ParameterLib::SpatialPosition const& pos,
        MatrixOfVelocityAtIntegrationPoints& darcy_velocity_at_ips) const;

    const int _gravitational_axis_id;
    const double _gravitational_acceleration;
    const double _reference_temperature;
    const LiquidFlowMaterialProperties& _material_properties;
};

}  // namespace LiquidFlow
}  // namespace ProcessLib

#include "LiquidFlowLocalAssembler-impl.h"
