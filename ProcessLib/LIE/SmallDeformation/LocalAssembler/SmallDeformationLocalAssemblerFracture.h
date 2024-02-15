/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "IntegrationPointDataFracture.h"
#include "NumLib/Fem/Integration/GenericIntegrationMethod.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/LIE/Common/FractureProperty.h"
#include "ProcessLib/LIE/Common/HMatrixUtils.h"
#include "ProcessLib/LIE/Common/JunctionProperty.h"
#include "ProcessLib/LIE/SmallDeformation/SmallDeformationProcessData.h"
#include "SecondaryData.h"
#include "SmallDeformationLocalAssemblerInterface.h"

namespace ProcessLib
{
namespace LIE
{
namespace SmallDeformation
{
template <typename ShapeFunction, int DisplacementDim>
class SmallDeformationLocalAssemblerFracture
    : public SmallDeformationLocalAssemblerInterface
{
public:
    using ShapeMatricesType =
        ShapeMatrixPolicyType<ShapeFunction, DisplacementDim>;
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;
    using HMatricesType = HMatrixPolicyType<ShapeFunction, DisplacementDim>;

    using HMatrixType = typename HMatricesType::HMatrixType;
    using StiffnessMatrixType = typename HMatricesType::StiffnessMatrixType;
    using NodalForceVectorType = typename HMatricesType::NodalForceVectorType;
    using NodalDisplacementVectorType =
        typename HMatricesType::NodalForceVectorType;

    using ForceVectorType = typename HMatricesType::ForceVectorType;
    using GlobalDimVectorType = Eigen::Matrix<double, DisplacementDim, 1>;

    SmallDeformationLocalAssemblerFracture(
        SmallDeformationLocalAssemblerFracture const&) = delete;
    SmallDeformationLocalAssemblerFracture(
        SmallDeformationLocalAssemblerFracture&&) = delete;

    SmallDeformationLocalAssemblerFracture(
        MeshLib::Element const& e,
        std::size_t const n_variables,
        std::size_t const local_matrix_size,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        SmallDeformationProcessData<DisplacementDim>& process_data);

    void assemble(double const /*t*/, double const /*dt*/,
                  std::vector<double> const& /*local_x*/,
                  std::vector<double> const& /*local_x_prev*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_b_data*/) override
    {
        OGS_FATAL(
            "SmallDeformationLocalAssembler: assembly without jacobian is not "
            "implemented.");
    }

    void assembleWithJacobian(double const t, double const dt,
                              Eigen::VectorXd const& local_u,
                              Eigen::VectorXd& local_b,
                              Eigen::MatrixXd& local_J) override;

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

    void computeSecondaryVariableConcreteWithVector(
        const double t, Eigen::VectorXd const& local_u) override;

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _secondary_data.N[integration_point];

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

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

private:
    SmallDeformationProcessData<DisplacementDim>& _process_data;
    std::vector<FractureProperty*> _fracture_props;
    std::vector<JunctionProperty*> _junction_props;
    std::unordered_map<int, int> _fracID_to_local;
    FractureProperty const* _fracture_property = nullptr;

    using IntegrationPointDataType =
        IntegrationPointDataFracture<HMatricesType, DisplacementDim>;
    std::vector<IntegrationPointDataType,
                Eigen::aligned_allocator<IntegrationPointDataType>>
        _ip_data;

    NumLib::GenericIntegrationMethod const& _integration_method;
    std::vector<ShapeMatrices, Eigen::aligned_allocator<
                                   typename ShapeMatricesType::ShapeMatrices>>
        _shape_matrices;
    MeshLib::Element const& _element;
    SecondaryData<typename ShapeMatrices::ShapeType> _secondary_data;
};

}  // namespace SmallDeformation
}  // namespace LIE
}  // namespace ProcessLib

#include "SmallDeformationLocalAssemblerFracture-impl.h"
