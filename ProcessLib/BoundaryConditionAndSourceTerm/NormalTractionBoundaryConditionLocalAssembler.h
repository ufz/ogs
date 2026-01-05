// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "GenericNaturalBoundaryConditionLocalAssembler.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "MeshLib/Elements/FaceRule.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
namespace NormalTractionBoundaryCondition
{
template <typename ShapeMatricesType>
struct IntegrationPointData final
{
    IntegrationPointData(
        typename ShapeMatricesType::ShapeMatrices::ShapeType const N_,
        typename ShapeMatricesType::GlobalDimVectorType const n_,
        double const integration_weight_)
        : N(N_), n(n_), integration_weight(integration_weight_)
    {
    }

    typename ShapeMatricesType::ShapeMatrices::ShapeType const N;
    typename ShapeMatricesType::GlobalDimVectorType const n;
    double const integration_weight;
};

class NormalTractionBoundaryConditionLocalAssemblerInterface
{
public:
    virtual void assemble(
        std::size_t const id,
        NumLib::LocalToGlobalIndexMap const& dof_table_boundary, double const t,
        std::vector<GlobalVector*> const& /*x*/, GlobalMatrix* /*K*/,
        GlobalVector& b, GlobalMatrix* /*Jac*/) = 0;
    virtual ~NormalTractionBoundaryConditionLocalAssemblerInterface() = default;
};

template <typename ShapeFunctionDisplacement, int GlobalDim>
class NormalTractionBoundaryConditionLocalAssembler final
    : public NormalTractionBoundaryConditionLocalAssemblerInterface
{
public:
    using ShapeMatricesType =
        ShapeMatrixPolicyType<ShapeFunctionDisplacement, GlobalDim>;
    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;

    NormalTractionBoundaryConditionLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const local_matrix_size,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        ParameterLib::Parameter<double> const& pressure,
        std::vector<Eigen::Vector3d> const& element_normals)
        : _integration_method(integration_method),
          _pressure(pressure),
          _element(e)
    {
        _local_rhs.setZero(local_matrix_size);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        _ip_data.reserve(n_integration_points);

        auto const shape_matrices_u =
            NumLib::initShapeMatrices<ShapeFunctionDisplacement,
                                      ShapeMatricesType, GlobalDim>(
                e, is_axially_symmetric, _integration_method);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            double const integration_weight =
                _integration_method.getWeightedPoint(ip).getWeight() *
                shape_matrices_u[ip].integralMeasure *
                shape_matrices_u[ip].detJ;

            _ip_data.emplace_back(shape_matrices_u[ip].N,
                                  element_normals[e.getID()].head<GlobalDim>(),
                                  integration_weight);
        }
    }

    void assemble(std::size_t const id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const t, std::vector<GlobalVector*> const& /*x*/,
                  GlobalMatrix* /*K*/, GlobalVector& local_rhs,
                  GlobalMatrix* /*Jac*/) override
    {
        _local_rhs.setZero();

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        NodalVectorType pressure =
            _pressure.getNodalValuesOnElement(_element, t);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto const& w = _ip_data[ip].integration_weight;
            auto const& N = _ip_data[ip].N;
            auto const& n = _ip_data[ip].n;

            typename ShapeMatricesType::template MatrixType<GlobalDim,
                                                            displacement_size>
                N_u = ShapeMatricesType::template MatrixType<
                    GlobalDim, displacement_size>::Zero(GlobalDim,
                                                        displacement_size);
            for (int i = 0; i < GlobalDim; ++i)
            {
                N_u.template block<1, displacement_size / GlobalDim>(
                       i, i * displacement_size / GlobalDim)
                    .noalias() = N;
            }

            _local_rhs.noalias() -= n.transpose() * N_u * pressure.dot(N) * w;
        }

        auto const indices = NumLib::getIndices(id, dof_table_boundary);
        local_rhs.add(indices, _local_rhs);
    }

private:
    NumLib::GenericIntegrationMethod const& _integration_method;
    ParameterLib::Parameter<double> const& _pressure;

    static const int displacement_size =
        ShapeFunctionDisplacement::NPOINTS * GlobalDim;
    std::vector<
        IntegrationPointData<ShapeMatricesType>,
        Eigen::aligned_allocator<IntegrationPointData<ShapeMatricesType>>>
        _ip_data;

    typename ShapeMatricesType::template VectorType<displacement_size>
        _local_rhs;

    MeshLib::Element const& _element;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace NormalTractionBoundaryCondition
}  // namespace ProcessLib
