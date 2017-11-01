/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <vector>

#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "LocalAssemblerInterface.h"
#include "SmallDeformationProcessData.h"

namespace ProcessLib
{
namespace SmallDeformation
{
template <typename BMatricesType, typename ShapeMatricesType,
          int DisplacementDim>
struct IntegrationPointData final
{
    explicit IntegrationPointData(
        MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material)
        : solid_material(solid_material),
          material_state_variables(
              solid_material.createMaterialStateVariables())
    {
    }

    typename BMatricesType::KelvinVectorType sigma, sigma_prev;
    typename BMatricesType::KelvinVectorType eps, eps_prev;

    MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        material_state_variables;

    double integration_weight;
    typename ShapeMatricesType::NodalRowVectorType N;
    typename ShapeMatricesType::GlobalDimNodalMatrixType dNdx;

    void pushBackState()
    {
        eps_prev = eps;
        sigma_prev = sigma;
        material_state_variables->pushBackState();
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

/// Used by for extrapolation of the integration point values. It is ordered
/// (and stored) by integration points.
template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType, Eigen::aligned_allocator<ShapeMatrixType>> N;
};

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
class SmallDeformationLocalAssembler
    : public SmallDeformationLocalAssemblerInterface<DisplacementDim>
{
public:
    using ShapeMatricesType =
        ShapeMatrixPolicyType<ShapeFunction, DisplacementDim>;
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;
    using BMatricesType = BMatrixPolicyType<ShapeFunction, DisplacementDim>;

    using BMatrixType = typename BMatricesType::BMatrixType;
    using StiffnessMatrixType = typename BMatricesType::StiffnessMatrixType;
    using NodalForceVectorType = typename BMatricesType::NodalForceVectorType;
    using NodalDisplacementVectorType =
        typename BMatricesType::NodalForceVectorType;

    SmallDeformationLocalAssembler(SmallDeformationLocalAssembler const&) =
        delete;
    SmallDeformationLocalAssembler(SmallDeformationLocalAssembler&&) = delete;

    SmallDeformationLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        SmallDeformationProcessData<DisplacementDim>& process_data)
        : _process_data(process_data),
          _integration_method(integration_order),
          _element(e),
          _is_axially_symmetric(is_axially_symmetric)
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        _ip_data.reserve(n_integration_points);
        _secondary_data.N.resize(n_integration_points);

        auto const shape_matrices =
            initShapeMatrices<ShapeFunction, ShapeMatricesType,
                              IntegrationMethod, DisplacementDim>(
                e, is_axially_symmetric, _integration_method);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data.emplace_back(*_process_data.material);
            auto& ip_data = _ip_data[ip];
            auto const& sm = shape_matrices[ip];
            _ip_data[ip].integration_weight =
                _integration_method.getWeightedPoint(ip).getWeight() *
                sm.integralMeasure * sm.detJ;

            ip_data.N = sm.N;
            ip_data.dNdx = sm.dNdx;

            // Initialize current time step values
            ip_data.sigma.setZero(
                KelvinVectorDimensions<DisplacementDim>::value);
            ip_data.eps.setZero(KelvinVectorDimensions<DisplacementDim>::value);

            // Previous time step values are not initialized and are set later.
            ip_data.sigma_prev.resize(
                KelvinVectorDimensions<DisplacementDim>::value);
            ip_data.eps_prev.resize(
                KelvinVectorDimensions<DisplacementDim>::value);

            _secondary_data.N[ip] = shape_matrices[ip].N;
        }
    }

    void assemble(double const /*t*/, std::vector<double> const& local_x,
            std::vector<double>& /*local_M_data*/,
            std::vector<double>& local_K_data,
            std::vector<double>& local_b_data) override
    {
        //only for linear elastic
        //        M=0;
        //        K = B.transpose() * C B * w //vormals Jac
        //        b = N_u_op.transpose() * rho * b * w;

        auto const local_matrix_size = local_x.size();

        auto local_K = MathLib::createZeroedMatrix<StiffnessMatrixType>(
            local_K_data, local_matrix_size, local_matrix_size);

        auto local_b = MathLib::createZeroedVector<NodalDisplacementVectorType>(
            local_b_data, local_matrix_size);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            KelvinMatrixType<DisplacementDim> C;
            std::tie(sigma, state, C) = std::move(*solution);

            auto const rho = _process_data.solid_density(t, x_position)[0];
            auto const& b = _process_data.specific_body_force;
            local_b.noalias() -=
                (B.transpose() * sigma - N_u_op.transpose() * rho * b) * w;
            local_K.noalias() += B.transpose() * C * B * w;
        }

    }

    void assembleWithJacobian(double const t,
                              std::vector<double> const& local_x,
                              std::vector<double> const& /*local_xdot*/,
                              const double /*dxdot_dx*/, const double /*dx_dx*/,
                              std::vector<double>& /*local_M_data*/,
                              std::vector<double>& /*local_K_data*/,
                              std::vector<double>& local_b_data,
                              std::vector<double>& local_Jac_data) override
    {
        auto const local_matrix_size = local_x.size();

        auto local_Jac = MathLib::createZeroedMatrix<StiffnessMatrixType>(
            local_Jac_data, local_matrix_size, local_matrix_size);

        auto local_b = MathLib::createZeroedVector<NodalDisplacementVectorType>(
            local_b_data, local_matrix_size);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto const& w = _ip_data[ip].integration_weight;
            auto const& N = _ip_data[ip].N;
            auto const& dNdx = _ip_data[ip].dNdx;

            typename ShapeMatricesType::template MatrixType<DisplacementDim,
                                                            displacement_size>
                N_u_op = ShapeMatricesType::template MatrixType<
                    DisplacementDim,
                    displacement_size>::Zero(DisplacementDim,
                                             displacement_size);
            for (int i = 0; i < DisplacementDim; ++i)
                N_u_op
                    .template block<1, displacement_size / DisplacementDim>(
                        i, i * displacement_size / DisplacementDim)
                    .noalias() = N;

            auto const x_coord =
                interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(
                    _element, N);
            auto const B = LinearBMatrix::computeBMatrix<
                DisplacementDim, ShapeFunction::NPOINTS,
                typename BMatricesType::BMatrixType>(dNdx, N, x_coord,
                                                     _is_axially_symmetric);

            auto const& eps_prev = _ip_data[ip].eps_prev;
            auto const& sigma_prev = _ip_data[ip].sigma_prev;

            auto& eps = _ip_data[ip].eps;
            auto& sigma = _ip_data[ip].sigma;
            auto& state = _ip_data[ip].material_state_variables;

            eps.noalias() =
                B *
                Eigen::Map<typename BMatricesType::NodalForceVectorType const>(
                    local_x.data(), ShapeFunction::NPOINTS * DisplacementDim);

            auto&& solution = _ip_data[ip].solid_material.integrateStress(
                t, x_position, _process_data.dt, eps_prev, eps, sigma_prev,
                *state);

            if (!solution)
                OGS_FATAL("Computation of local constitutive relation failed.");

            KelvinMatrixType<DisplacementDim> C;
            std::tie(sigma, state, C) = std::move(*solution);

            auto const rho = _process_data.solid_density(t, x_position)[0];
            auto const& b = _process_data.specific_body_force;
            local_b.noalias() -=
                (B.transpose() * sigma - N_u_op.transpose() * rho * b) * w;
            local_Jac.noalias() += B.transpose() * C * B * w;
        }
    }

    void preTimestepConcrete(std::vector<double> const& /*local_x*/,
                             double const /*t*/,
                             double const /*delta_t*/) override
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

    std::vector<double> const& getIntPtSigma(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        using KelvinVectorType = typename BMatricesType::KelvinVectorType;
        auto const kelvin_vector_size =
            KelvinVectorDimensions<DisplacementDim>::value;
        auto const num_intpts = _ip_data.size();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, kelvin_vector_size, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, kelvin_vector_size, num_intpts);

        // TODO make a general implementation for converting KelvinVectors
        // back to symmetric rank-2 tensors.
        for (unsigned ip = 0; ip < num_intpts; ++ip)
        {
            auto const& sigma = _ip_data[ip].sigma;

            for (typename KelvinVectorType::Index component = 0;
                 component < kelvin_vector_size && component < 3;
                 ++component)
            {  // xx, yy, zz components
                cache_mat(component, ip) = sigma[component];
            }
            for (typename KelvinVectorType::Index component = 3;
                 component < kelvin_vector_size;
                 ++component)
            {  // mixed xy, yz, xz components
                cache_mat(component, ip) = sigma[component] / std::sqrt(2);
            }
        }

        return cache;
    }

    // TODO remove the component-wise methods.
    std::vector<double> const& getIntPtSigmaXX(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtSigma(cache, 0);
    }

    std::vector<double> const& getIntPtSigmaYY(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtSigma(cache, 1);
    }

    std::vector<double> const& getIntPtSigmaZZ(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtSigma(cache, 2);
    }

    std::vector<double> const& getIntPtSigmaXY(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtSigma(cache, 3);
    }

    std::vector<double> const& getIntPtSigmaYZ(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtSigma(cache, 4);
    }

    std::vector<double> const& getIntPtSigmaXZ(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtSigma(cache, 5);
    }

    virtual std::vector<double> const& getIntPtEpsilon(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        using KelvinVectorType = typename BMatricesType::KelvinVectorType;
        auto const kelvin_vector_size =
            KelvinVectorDimensions<DisplacementDim>::value;
        auto const num_intpts = _ip_data.size();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, kelvin_vector_size, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, kelvin_vector_size, num_intpts);

        // TODO make a general implementation for converting KelvinVectors
        // back to symmetric rank-2 tensors.
        for (unsigned ip = 0; ip < num_intpts; ++ip)
        {
            auto const& eps = _ip_data[ip].eps;

            for (typename KelvinVectorType::Index component = 0;
                 component < kelvin_vector_size && component < 3;
                 ++component)
            {  // xx, yy, zz components
                cache_mat(component, ip) = eps[component];
            }
            for (typename KelvinVectorType::Index component = 3;
                 component < kelvin_vector_size;
                 ++component)
            {  // mixed xy, yz, xz components
                cache_mat(component, ip) = eps[component] / std::sqrt(2);
            }
        }

        return cache;
    }

    // TODO remove the component-wise methods
    std::vector<double> const& getIntPtEpsilonXX(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtEpsilon(cache, 0);
    }

    std::vector<double> const& getIntPtEpsilonYY(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtEpsilon(cache, 1);
    }

    std::vector<double> const& getIntPtEpsilonZZ(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtEpsilon(cache, 2);
    }

    std::vector<double> const& getIntPtEpsilonXY(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtEpsilon(cache, 3);
    }

    std::vector<double> const& getIntPtEpsilonYZ(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtEpsilon(cache, 4);
    }

    std::vector<double> const& getIntPtEpsilonXZ(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtEpsilon(cache, 5);
    }

    unsigned getNumberOfIntegrationPoints() const override
    {
        return _integration_method.getNumberOfPoints();
    }

    typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables const&
    getMaterialStateVariablesAt(unsigned integration_point) const override
    {
        return *_ip_data[integration_point].material_state_variables;
    }

private:
    std::vector<double> const& getIntPtSigma(std::vector<double>& cache,
                                             std::size_t const component) const
    {
        cache.clear();
        cache.reserve(_ip_data.size());

        for (auto const& ip_data : _ip_data)
        {
            if (component < 3)  // xx, yy, zz components
                cache.push_back(ip_data.sigma[component]);
            else  // mixed xy, yz, xz components
                cache.push_back(ip_data.sigma[component] / std::sqrt(2));
        }

        return cache;
    }

    std::vector<double> const& getIntPtEpsilon(
        std::vector<double>& cache, std::size_t const component) const
    {
        cache.clear();
        cache.reserve(_ip_data.size());

        for (auto const& ip_data : _ip_data)
        {
            if (component < 3)  // xx, yy, zz components
                cache.push_back(ip_data.eps[component]);
            else  // mixed xy, yz, xz components
                cache.push_back(ip_data.eps[component] / std::sqrt(2));
        }

        return cache;
    }

    SmallDeformationProcessData<DisplacementDim>& _process_data;

    std::vector<
        IntegrationPointData<BMatricesType, ShapeMatricesType, DisplacementDim>,
        Eigen::aligned_allocator<IntegrationPointData<
            BMatricesType, ShapeMatricesType, DisplacementDim>>>
        _ip_data;

    IntegrationMethod _integration_method;
    MeshLib::Element const& _element;
    SecondaryData<typename ShapeMatrices::ShapeType> _secondary_data;
    bool const _is_axially_symmetric;

    static const int displacement_size =
        ShapeFunction::NPOINTS * DisplacementDim;
};

}  // namespace SmallDeformation
}  // namespace ProcessLib
