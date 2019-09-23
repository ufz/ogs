/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file
 *  Created on October 11, 2017, 1:42 PM
 */

#pragma once

#include <Eigen/Dense>
#include <vector>

#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "HTFEM.h"

namespace ProcessLib
{
struct LocalCoupledSolutions;

namespace HT
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class StaggeredHTFEM : public HTFEM<ShapeFunction, IntegrationMethod, GlobalDim>
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalMatrixType =
        typename ShapeMatricesType::template MatrixType<ShapeFunction::NPOINTS,
                                                        ShapeFunction::NPOINTS>;
    using LocalVectorType =
        typename ShapeMatricesType::template VectorType<ShapeFunction::NPOINTS>;

    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using NodalRowVectorType = typename ShapeMatricesType::NodalRowVectorType;

    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;
    using GlobalDimNodalMatrixType =
        typename ShapeMatricesType::GlobalDimNodalMatrixType;
    using GlobalDimMatrixType = typename ShapeMatricesType::GlobalDimMatrixType;

public:
    StaggeredHTFEM(MeshLib::Element const& element,
                   std::size_t const local_matrix_size,
                   bool is_axially_symmetric,
                   unsigned const integration_order,
                   HTProcessData const& process_data,
                   const int heat_transport_process_id,
                   const int hydraulic_process_id)
        : HTFEM<ShapeFunction, IntegrationMethod, GlobalDim>(
              element, local_matrix_size, is_axially_symmetric,
              integration_order, process_data, 1),
        _heat_transport_process_id(heat_transport_process_id),
        _hydraulic_process_id(hydraulic_process_id)
    {
    }

    void assembleForStaggeredScheme(
        double const t, std::vector<double>& local_M_data,
        std::vector<double>& local_K_data, std::vector<double>& local_b_data,
        LocalCoupledSolutions const& coupled_xs) override;

    std::vector<double> const& getIntPtDarcyVelocity(
        const double t,
        GlobalVector const& current_solution,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        std::vector<double>& cache) const override;

private:
    void assembleHydraulicEquation(double const t,
                                   std::vector<double>& local_M_data,
                                   std::vector<double>& local_K_data,
                                   std::vector<double>& local_b_data,
                                   LocalCoupledSolutions const& coupled_xs);

    void assembleHeatTransportEquation(
        double const t, std::vector<double>& local_M_data,
        std::vector<double>& local_K_data, std::vector<double>& local_b_data,
        LocalCoupledSolutions const& coupled_xs);
    const int _heat_transport_process_id;
    const int _hydraulic_process_id;
};

}  // namespace HT
}  // namespace ProcessLib

#include "StaggeredHTFEM-impl.h"
