/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "MeshLib/Elements/FaceRule.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "ParameterLib/Parameter.h"

#include "GenericNaturalBoundaryConditionLocalAssembler.h"

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
        std::vector<GlobalVector*> const& /*x*/, GlobalMatrix& /*K*/,
        GlobalVector& b, GlobalMatrix* /*Jac*/) = 0;
    virtual ~NormalTractionBoundaryConditionLocalAssemblerInterface() = default;
};

template <typename ShapeFunctionDisplacement, typename IntegrationMethod,
          unsigned GlobalDim>
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
        bool const is_axially_symmetric,
        unsigned const integration_order,
        ParameterLib::Parameter<double> const& pressure)
        : integration_method_(integration_order),
          pressure_(pressure),
          element_(e)
    {
        local_rhs_.setZero(local_matrix_size);

        unsigned const n_integration_points =
            integration_method_.getNumberOfPoints();

        ip_data_.reserve(n_integration_points);

        auto const shape_matrices_u =
            initShapeMatrices<ShapeFunctionDisplacement, ShapeMatricesType,
                              IntegrationMethod, GlobalDim>(
                e, is_axially_symmetric, integration_method_);

        GlobalDimVectorType element_normal(GlobalDim);

        // TODO Extend to rotated 2d meshes and line elements.
        if (e.getGeomType() == MeshLib::MeshElemType::LINE)
        {
            auto v1 = (*e.getNode(1)) - (*e.getNode(0));
            element_normal[0] = -v1[1];
            element_normal[1] = v1[0];
            element_normal.normalize();
        }
        else
        {
            auto const element_normal_vector =
                MeshLib::FaceRule::getSurfaceNormal(&e).getNormalizedVector();

            std::copy_n(element_normal_vector.getCoords(), GlobalDim,
                        element_normal.data());
        }

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {

            double const integration_weight =
                integration_method_.getWeightedPoint(ip).getWeight() *
                shape_matrices_u[ip].integralMeasure *
                shape_matrices_u[ip].detJ;

            ip_data_.emplace_back(shape_matrices_u[ip].N, element_normal,
                                  integration_weight);
        }
    }

    void assemble(std::size_t const id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const t, std::vector<GlobalVector*> const& /*x*/,
                  GlobalMatrix& /*K*/, GlobalVector& local_rhs,
                  GlobalMatrix* /*Jac*/) override
    {
        local_rhs_.setZero();

        unsigned const n_integration_points =
            integration_method_.getNumberOfPoints();

        NodalVectorType pressure =
            pressure_.getNodalValuesOnElement(element_, t);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto const& w = ip_data_[ip].integration_weight;
            auto const& N = ip_data_[ip].N;
            auto const& n = ip_data_[ip].n;

            typename ShapeMatricesType::template MatrixType<GlobalDim,
                                                            displacement_size>
                N_u = ShapeMatricesType::template MatrixType<
                    GlobalDim, displacement_size>::Zero(GlobalDim,
                                                        displacement_size);
            for (int i = 0; i < static_cast<int>(GlobalDim); ++i)
            {
                N_u.template block<1, displacement_size / GlobalDim>(
                       i, i * displacement_size / GlobalDim)
                    .noalias() = N;
            }

            local_rhs_.noalias() -= n.transpose() * N_u * pressure.dot(N) * w;
        }

        auto const indices = NumLib::getIndices(id, dof_table_boundary);
        local_rhs.add(indices, local_rhs_);
    }

private:
    IntegrationMethod const integration_method_;
    ParameterLib::Parameter<double> const& pressure_;

    static const int displacement_size =
        ShapeFunctionDisplacement::NPOINTS * GlobalDim;
    std::vector<
        IntegrationPointData<ShapeMatricesType>,
        Eigen::aligned_allocator<IntegrationPointData<ShapeMatricesType>>>
        ip_data_;

    typename ShapeMatricesType::template VectorType<displacement_size>
        local_rhs_;

    MeshLib::Element const& element_;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace NormalTractionBoundaryCondition
}  // namespace ProcessLib
