/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * Created on August 19, 2016, 2:28 PM
 */

#pragma once

#include <Eigen/Core>
#include <vector>

#include "LiquidFlowData.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEffectiveThermalConductivity.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/Integration/GenericIntegrationMethod.h"
#include "NumLib/Fem/ShapeMatrixCache.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"

namespace ProcessLib
{
namespace LiquidFlow
{
template <typename NodalRowVectorType, typename GlobalDimNodalMatrixType>
struct IntegrationPointData final
{
    explicit IntegrationPointData(GlobalDimNodalMatrixType const& dNdx_,
                                  double const& integration_weight_)
        : dNdx(dNdx_), integration_weight(integration_weight_)
    {
    }

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
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_tables,
        std::vector<double>& cache) const = 0;
};

template <typename ShapeFunction, int GlobalDim>
class LiquidFlowLocalAssembler : public LiquidFlowLocalAssemblerInterface
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalAssemblerTraits = ProcessLib::LocalAssemblerTraits<
        ShapeMatricesType, ShapeFunction::NPOINTS, NUM_NODAL_DOF, GlobalDim>;

    using NodalMatrixType = typename LocalAssemblerTraits::LocalMatrix;
    using NodalVectorType = typename LocalAssemblerTraits::LocalVector;
    using NodalRowVectorType = typename ShapeMatricesType::NodalRowVectorType;
    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;
    using GlobalDimMatrixType = typename ShapeMatricesType::GlobalDimMatrixType;
    using GlobalDimNodalMatrixType =
        typename ShapeMatricesType::GlobalDimNodalMatrixType;

public:
    LiquidFlowLocalAssembler(
        MeshLib::Element const& element,
        std::size_t const /*local_matrix_size*/,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        LiquidFlowData const& process_data,
        NumLib::ShapeMatrixCache const& shape_matrix_cache)
        : _element(element),
          _integration_method(integration_method),
          _process_data(process_data),
          _shape_matrix_cache(shape_matrix_cache)
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        _ip_data.reserve(n_integration_points);

        auto const& shape_matrices =
            NumLib::initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                      GlobalDim>(element, is_axially_symmetric,
                                                 _integration_method);

        ParameterLib::SpatialPosition pos;
        pos.setElementID(_element.getID());

        double const aperture_size =
            (_element.getDimension() == 3u)
                ? 1.0
                : _process_data.aperture_size(0.0, pos)[0];

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data.emplace_back(
                shape_matrices[ip].dNdx,
                _integration_method.getWeightedPoint(ip).getWeight() *
                    shape_matrices[ip].integralMeasure *
                    shape_matrices[ip].detJ * aperture_size);
        }
    }

    void assemble(double const t, double const dt,
                  std::vector<double> const& local_x,
                  std::vector<double> const& /*local_x_prev*/,
                  std::vector<double>& local_M_data,
                  std::vector<double>& local_K_data,
                  std::vector<double>& local_b_data) override;

    /// Computes the flux in the point \c p_local_coords that is given in local
    /// coordinates using the values from \c local_x.
    Eigen::Vector3d getFlux(MathLib::Point3d const& p_local_coords,
                            double const t,
                            std::vector<double> const& local_x) const override;

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N =
            _shape_matrix_cache
                .NsHigherOrder<typename ShapeFunction::MeshElement>();

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(
            N[integration_point].data(), N[integration_point].size());
    }

    std::vector<double> const& getIntPtDarcyVelocity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& velocity_cache) const override;

private:
    MeshLib::Element const& _element;

    NumLib::GenericIntegrationMethod const& _integration_method;
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
            GlobalDimMatrixType const& permeability, double const mu,
            double const rho_L, GlobalDimVectorType const& specific_body_force,
            bool const has_gravity);

        static Eigen::Matrix<double, GlobalDim, 1> calculateVelocity(
            Eigen::Map<const NodalVectorType> const& local_p,
            IntegrationPointData<NodalRowVectorType,
                                 GlobalDimNodalMatrixType> const& ip_data,
            GlobalDimMatrixType const& permeability, double const mu,
            double const rho_L, GlobalDimVectorType const& specific_body_force,
            bool const has_gravity);
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
            GlobalDimMatrixType const& permeability, double const mu,
            double const rho_L, GlobalDimVectorType const& specific_body_force,
            bool const has_gravity);

        static Eigen::Matrix<double, GlobalDim, 1> calculateVelocity(
            Eigen::Map<const NodalVectorType> const& local_p,
            IntegrationPointData<NodalRowVectorType,
                                 GlobalDimNodalMatrixType> const& ip_data,
            GlobalDimMatrixType const& permeability, double const mu,
            double const rho_L, GlobalDimVectorType const& specific_body_force,
            bool const has_gravity);
    };

    template <typename LaplacianGravityVelocityCalculator>
    void assembleMatrixAndVector(double const t, double const dt,
                                 std::vector<double> const& local_x,
                                 std::vector<double>& local_M_data,
                                 std::vector<double>& local_K_data,
                                 std::vector<double>& local_b_data);

    template <typename LaplacianGravityVelocityCalculator,
              typename VelocityCacheType>
    void computeProjectedDarcyVelocity(
        const double t, const double dt, std::vector<double> const& local_x,
        ParameterLib::SpatialPosition const& pos,
        VelocityCacheType& darcy_velocity_at_ips) const;

    template <typename VelocityCacheType>
    void computeDarcyVelocity(bool const is_scalar_permeability, const double t,
                              const double dt,
                              std::vector<double> const& local_x,
                              ParameterLib::SpatialPosition const& pos,
                              VelocityCacheType& darcy_velocity_at_ips) const;

    const LiquidFlowData& _process_data;
    NumLib::ShapeMatrixCache const& _shape_matrix_cache;
};

}  // namespace LiquidFlow
}  // namespace ProcessLib

#include "LiquidFlowLocalAssembler-impl.h"
