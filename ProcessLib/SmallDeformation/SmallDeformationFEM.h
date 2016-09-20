/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_SMALLDEFORMATION_FEM_H_
#define PROCESS_LIB_SMALLDEFORMATION_FEM_H_

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

#include "SmallDeformationProcessData.h"

template <typename BMatricesType, int DisplacementDim>
struct IntegrationPointData final
{
    explicit IntegrationPointData(
        Solids::MechanicsBase<DisplacementDim>& solid_material)
        : _solid_material(solid_material),
          _material_state_variables(
              _solid_material.createMaterialStateVariables())
    {
    }

#if defined(_MSC_VER) && _MSC_VER < 1900
    // The default generated move-ctor is correctly generated for other
    // compilers.
    explicit IntegrationPointData(IntegrationPointData&& other)
        : _b_matrices(std::move(other._b_matrices)),
          _sigma(std::move(other._sigma)),
          _sigma_prev(std::move(other._sigma_prev)),
          _eps(std::move(other._eps)),
          _eps_prev(std::move(other._eps_prev)),
          _solid_material(other._solid_material),
          _material_state_variables(std::move(other._material_state_variables)),
          _C(std::move(other._C)),
          _detJ(std::move(other._detJ))
    {
    }
#endif  // _MSC_VER

    typename BMatricesType::BMatrixType _b_matrices;
    typename BMatricesType::KelvinVectorType _sigma, _sigma_prev;
    typename BMatricesType::KelvinVectorType _eps, _eps_prev;

    Solids::MechanicsBase<DisplacementDim>& _solid_material;
    std::unique_ptr<
        typename Solids::MechanicsBase<DisplacementDim>::MaterialStateVariables>
        _material_state_variables;

    typename BMatricesType::KelvinMatrixType _C;
    double _detJ;

    void pushBackState()
    {
        _eps_prev = _eps;
        _sigma_prev = _sigma;
        _material_state_variables->pushBackState();
    }
};

namespace ProcessLib
{
namespace SmallDeformation
{

/// Used by for extrapolation of the integration point values. It is ordered
/// (and stored) by integration points.
template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType> N;
    std::vector<double> _sigmaXX;
    std::vector<double> _sigmaYY;
    std::vector<double> _sigmaXY;
};

struct SmallDeformationLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
    virtual std::vector<double> const& getIntPtSigmaXX(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtSigmaYY(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtSigmaXY(
        std::vector<double>& /*cache*/) const = 0;
};

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
class SmallDeformationLocalAssembler
    : public SmallDeformationLocalAssemblerInterface
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
        unsigned const integration_order,
        SmallDeformationProcessData<DisplacementDim>& process_data)
        : _process_data(process_data),
          _integration_method(integration_order),
          _element(e)
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        _ip_data.reserve(n_integration_points);

        _secondary_data.N.resize(n_integration_points);
        _secondary_data._sigmaXX.resize(n_integration_points);
        _secondary_data._sigmaYY.resize(n_integration_points);
        _secondary_data._sigmaXY.resize(n_integration_points);

        auto const shape_matrices =
            initShapeMatrices<ShapeFunction, ShapeMatricesType,
                              IntegrationMethod, DisplacementDim>(
                e, _integration_method);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data.emplace_back(*_process_data._material);
            auto& ip_data = _ip_data[ip];
            ip_data._detJ = shape_matrices[ip].detJ;
            ip_data._b_matrices.resize(
                KelvinVectorDimensions<DisplacementDim>::value,
                ShapeFunction::NPOINTS * DisplacementDim);
            LinearBMatrix::computeBMatrix<
                DisplacementDim, ShapeFunction::NPOINTS,
                typename ShapeMatricesType::GlobalDimNodalMatrixType,
                BMatrixType>(shape_matrices[ip].dNdx, ip_data._b_matrices);

            ip_data._sigma.resize(
                KelvinVectorDimensions<DisplacementDim>::value);
            ip_data._sigma_prev.resize(
                KelvinVectorDimensions<DisplacementDim>::value);
            ip_data._eps.resize(KelvinVectorDimensions<DisplacementDim>::value);
            ip_data._eps_prev.resize(
                KelvinVectorDimensions<DisplacementDim>::value);
            ip_data._C.resize(KelvinVectorDimensions<DisplacementDim>::value,
                              KelvinVectorDimensions<DisplacementDim>::value);

            _secondary_data.N[ip] = shape_matrices[ip].N;
        }
    }

    void assemble(double const /*t*/, std::vector<double> const& /*local_x*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_b_data*/) override
    {
        OGS_FATAL(
            "SmallDeformationLocalAssembler: assembly without jacobian is not "
            "implemented.");
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
            auto const& wp = _integration_method.getWeightedPoint(ip);
            auto const& detJ = _ip_data[ip]._detJ;

            auto const& B = _ip_data[ip]._b_matrices;
            auto const& eps_prev = _ip_data[ip]._eps_prev;
            auto const& sigma_prev = _ip_data[ip]._sigma_prev;

            auto& eps = _ip_data[ip]._eps;
            auto& sigma = _ip_data[ip]._sigma;
            auto& C = _ip_data[ip]._C;
            auto& material_state_variables =
                *_ip_data[ip]._material_state_variables;

            eps.noalias() =
                B *
                Eigen::Map<typename BMatricesType::NodalForceVectorType const>(
                    local_x.data(), ShapeFunction::NPOINTS * DisplacementDim);

            _ip_data[ip]._solid_material.computeConstitutiveRelation(
                t, x_position, _process_data.dt, eps_prev, eps, sigma_prev,
                sigma, C, material_state_variables);

            local_b.noalias() -= B.transpose() * sigma * detJ * wp.getWeight();
            local_Jac.noalias() +=
                B.transpose() * C * B * detJ * wp.getWeight();


            // TODO: Reuse _ip_data[ip]._sigma
            _secondary_data._sigmaXX[ip] = sigma[0];
            _secondary_data._sigmaYY[ip] = sigma[1];
            _secondary_data._sigmaXY[ip] = sigma[3];
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

    std::vector<double> const& getIntPtSigmaXX(
        std::vector<double>& /*cache*/) const override
    {
        return _secondary_data._sigmaXX;
    }

    std::vector<double> const& getIntPtSigmaYY(
        std::vector<double>& /*cache*/) const override
    {
        return _secondary_data._sigmaYY;
    }

    std::vector<double> const& getIntPtSigmaXY(
        std::vector<double>& /*cache*/) const override
    {
        return _secondary_data._sigmaXY;
    }

private:
    SmallDeformationProcessData<DisplacementDim>& _process_data;

    std::vector<IntegrationPointData<BMatricesType, DisplacementDim>> _ip_data;

    IntegrationMethod _integration_method;
    MeshLib::Element const& _element;
    SecondaryData<typename ShapeMatrices::ShapeType> _secondary_data;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim, int DisplacementDim>
class LocalAssemblerData final
    : public SmallDeformationLocalAssembler<ShapeFunction, IntegrationMethod,
                                            DisplacementDim>
{
public:
    LocalAssemblerData(LocalAssemblerData const&) = delete;
    LocalAssemblerData(LocalAssemblerData&&) = delete;

    LocalAssemblerData(
        MeshLib::Element const& e,
        std::size_t const local_matrix_size,
        unsigned const integration_order,
        SmallDeformationProcessData<DisplacementDim>& process_data)
        : SmallDeformationLocalAssembler<ShapeFunction, IntegrationMethod,
                                         DisplacementDim>(
              e, local_matrix_size, integration_order, process_data)
    {
    }
};

}  // namespace SmallDeformation
}  // namespace ProcessLib

#endif  // PROCESS_LIB_SMALLDEFORMATION_FEM_H_
