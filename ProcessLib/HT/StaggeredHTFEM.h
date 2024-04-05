/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  Created on October 11, 2017, 1:42 PM
 */

#pragma once

#include <Eigen/Core>
#include <vector>

#include "HTFEM.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/Parameter.h"

namespace ProcessLib
{

namespace HT
{
template <typename ShapeFunction, int GlobalDim>
class StaggeredHTFEM : public HTFEM<ShapeFunction, GlobalDim>
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalMatrixType =
        typename ShapeMatricesType::template MatrixType<ShapeFunction::NPOINTS,
                                                        ShapeFunction::NPOINTS>;
    using LocalVectorType =
        typename ShapeMatricesType::template VectorType<ShapeFunction::NPOINTS>;

    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalRowVectorType = typename ShapeMatricesType::NodalRowVectorType;

    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;
    using GlobalDimNodalMatrixType =
        typename ShapeMatricesType::GlobalDimNodalMatrixType;
    using GlobalDimMatrixType = typename ShapeMatricesType::GlobalDimMatrixType;

    using HTFEM<ShapeFunction, GlobalDim>::pressure_index;
    using HTFEM<ShapeFunction, GlobalDim>::pressure_size;
    using HTFEM<ShapeFunction, GlobalDim>::temperature_index;
    using HTFEM<ShapeFunction, GlobalDim>::temperature_size;

public:
    StaggeredHTFEM(MeshLib::Element const& element,
                   std::size_t const local_matrix_size,
                   NumLib::GenericIntegrationMethod const& integration_method,
                   bool is_axially_symmetric,
                   HTProcessData const& process_data)
        : HTFEM<ShapeFunction, GlobalDim>(element, local_matrix_size,
                                          integration_method,
                                          is_axially_symmetric, process_data, 1)
    {
    }

    void assembleForStaggeredScheme(double const t, double const dt,
                                    Eigen::VectorXd const& local_x,
                                    Eigen::VectorXd const& local_x_prev,
                                    int const process_id,
                                    std::vector<double>& local_M_data,
                                    std::vector<double>& local_K_data,
                                    std::vector<double>& local_b_data) override;

    std::vector<double> const& getIntPtDarcyVelocity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override;

private:
    void assembleHydraulicEquation(double const t, double const dt,
                                   Eigen::VectorXd const& local_x,
                                   Eigen::VectorXd const& local_x_prev,
                                   std::vector<double>& local_M_data,
                                   std::vector<double>& local_K_data,
                                   std::vector<double>& local_b_data);

    void assembleHeatTransportEquation(double const t,
                                       double const dt,
                                       Eigen::VectorXd const& local_x,
                                       std::vector<double>& local_M_data,
                                       std::vector<double>& local_K_data);
};

}  // namespace HT
}  // namespace ProcessLib

#include "StaggeredHTFEM-impl.h"
