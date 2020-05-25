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

#include "ProcessLib/HeatTransportBHE/HeatTransportBHEProcessData.h"

#include "HeatTransportBHEProcessAssemblerInterface.h"
#include "IntegrationPointDataBHE.h"
#include "SecondaryData.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
template <typename ShapeFunction, typename IntegrationMethod, typename BHEType>
class HeatTransportBHELocalAssemblerBHE
    : public HeatTransportBHELocalAssemblerInterface
{
    static constexpr int bhe_unknowns = BHEType::number_of_unknowns;
    static constexpr int single_bhe_unknowns_size = ShapeFunction::NPOINTS;
    static constexpr int soil_temperature_size = ShapeFunction::NPOINTS;
    static constexpr int soil_temperature_index = 0;
    static constexpr int bhe_unknowns_size =
        single_bhe_unknowns_size * bhe_unknowns;
    static constexpr int bhe_unknowns_index = ShapeFunction::NPOINTS;
    static constexpr int local_matrix_size =
        soil_temperature_size + bhe_unknowns_size;

public:
    using ShapeMatricesType =
        ShapeMatrixPolicyType<ShapeFunction, 3 /* GlobalDim */>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;
    // Using dynamic size, because the number of unknowns in the BHE is runtime
    // value.
    using BheLocalMatrixType =
        typename ShapeMatricesType::template MatrixType<local_matrix_size,
                                                        local_matrix_size>;
    HeatTransportBHELocalAssemblerBHE(
        HeatTransportBHELocalAssemblerBHE const&) = delete;
    HeatTransportBHELocalAssemblerBHE(HeatTransportBHELocalAssemblerBHE&&) =
        delete;

    HeatTransportBHELocalAssemblerBHE(
        MeshLib::Element const& e,
        BHEType const& bhe,
        bool const is_axially_symmetric,
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
        IntegrationPointDataBHE<ShapeMatricesType>,
        Eigen::aligned_allocator<IntegrationPointDataBHE<ShapeMatricesType>>>
        ip_data_;

    IntegrationMethod integration_method_;

    BHEType const& bhe_;

    std::size_t const element_id_;

    SecondaryData<typename ShapeMatrices::ShapeType> secondary_data_;

    Eigen::Vector3d element_direction_;

    typename ShapeMatricesType::template MatrixType<bhe_unknowns_size,
                                                    bhe_unknowns_size>
        R_matrix_;

    typename ShapeMatricesType::template MatrixType<soil_temperature_size,
                                                    soil_temperature_size>
        R_s_matrix_;

    typename ShapeMatricesType::template MatrixType<bhe_unknowns_size,
                                                    soil_temperature_size>
        R_pi_s_matrix_;
};
}  // namespace HeatTransportBHE
}  // namespace ProcessLib

#include "HeatTransportBHELocalAssemblerBHE-impl.h"
