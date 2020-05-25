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

#include <vector>

#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "ProcessLib/HeatTransportBHE/HeatTransportBHEProcessData.h"

#include "HeatTransportBHEProcessAssemblerInterface.h"
#include "IntegrationPointDataSoil.h"
#include "SecondaryData.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
template <typename ShapeFunction, typename IntegrationMethod>
class HeatTransportBHELocalAssemblerSoil
    : public HeatTransportBHELocalAssemblerInterface
{
public:
    using ShapeMatricesType =
        ShapeMatrixPolicyType<ShapeFunction, 3 /* GlobalDim */>;
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using NodalRowVectorType = typename ShapeMatricesType::NodalRowVectorType;

    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;
    using GlobalDimNodalMatrixType =
        typename ShapeMatricesType::GlobalDimNodalMatrixType;

    HeatTransportBHELocalAssemblerSoil(
        HeatTransportBHELocalAssemblerSoil const&) = delete;
    HeatTransportBHELocalAssemblerSoil(HeatTransportBHELocalAssemblerSoil&&) =
        delete;

    HeatTransportBHELocalAssemblerSoil(
        MeshLib::Element const& e,
        bool is_axially_symmetric,
        unsigned const integration_order,
        HeatTransportBHEProcessData& process_data);

    void assemble(double const /*t*/, double const /*dt*/,
                  std::vector<double> const& /*local_x*/,
                  std::vector<double> const& /*local_xdot*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_b_data*/) override;

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = secondary_data_.N[integration_point];

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

private:
    HeatTransportBHEProcessData& process_data_;

    std::vector<
        IntegrationPointDataSoil<NodalRowVectorType, GlobalDimNodalMatrixType>,
        Eigen::aligned_allocator<IntegrationPointDataSoil<
            NodalRowVectorType, GlobalDimNodalMatrixType>>>
        ip_data_;

    IntegrationMethod const integration_method_;

    std::vector<ShapeMatrices, Eigen::aligned_allocator<ShapeMatrices>>
        shape_matrices_;

    std::size_t const element_id_;

    SecondaryData<typename ShapeMatrices::ShapeType> secondary_data_;
};
}  // namespace HeatTransportBHE
}  // namespace ProcessLib

#include "HeatTransportBHELocalAssemblerSoil-impl.h"
