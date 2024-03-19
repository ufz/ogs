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

#include "HeatTransportBHEProcessAssemblerInterface.h"
#include "IntegrationPointDataBHE.h"
#include "NumLib/Fem/Integration/GenericIntegrationMethod.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/HeatTransportBHE/HeatTransportBHEProcessData.h"
#include "SecondaryData.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
template <typename ShapeFunction, typename BHEType>
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
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    // Using dynamic size, because the number of unknowns in the BHE is runtime
    // value.
    using BheLocalMatrixType =
        typename ShapeMatricesType::template MatrixType<local_matrix_size,
                                                        local_matrix_size>;
    using BheLocalVectorType =
        typename ShapeMatricesType::template VectorType<local_matrix_size>;

    HeatTransportBHELocalAssemblerBHE(
        HeatTransportBHELocalAssemblerBHE const&) = delete;
    HeatTransportBHELocalAssemblerBHE(HeatTransportBHELocalAssemblerBHE&&) =
        delete;

    HeatTransportBHELocalAssemblerBHE(
        MeshLib::Element const& e,
        NumLib::GenericIntegrationMethod const& integration_method,
        BHEType const& bhe,
        bool const is_axially_symmetric,
        HeatTransportBHEProcessData& process_data);

    void assemble(double const /*t*/, double const /*dt*/,
                  std::vector<double> const& /*local_x*/,
                  std::vector<double> const& /*local_x_prev*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_b_data*/) override;

    void assembleWithJacobian(double const t, double const dt,
                              std::vector<double> const& local_x,
                              std::vector<double> const& local_x_prev,
                              std::vector<double>& /*local_M_data*/,
                              std::vector<double>& /*local_K_data*/,
                              std::vector<double>& local_rhs_data,
                              std::vector<double>& local_Jac_data) override;

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _secondary_data.N[integration_point];

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

private:
    HeatTransportBHEProcessData& _process_data;

    std::vector<
        IntegrationPointDataBHE<ShapeMatricesType>,
        Eigen::aligned_allocator<IntegrationPointDataBHE<ShapeMatricesType>>>
        _ip_data;

    NumLib::GenericIntegrationMethod const& _integration_method;

    BHEType const& _bhe;

    std::size_t const _element_id;

    SecondaryData<typename ShapeMatrices::ShapeType> _secondary_data;

    Eigen::Vector3d _element_direction;

    typename ShapeMatricesType::template MatrixType<bhe_unknowns_size,
                                                    bhe_unknowns_size>
        _R_matrix;

    typename ShapeMatricesType::template MatrixType<soil_temperature_size,
                                                    soil_temperature_size>
        _R_s_matrix;

    typename ShapeMatricesType::template MatrixType<bhe_unknowns_size,
                                                    soil_temperature_size>
        _R_pi_s_matrix;
};
}  // namespace HeatTransportBHE
}  // namespace ProcessLib

#include "HeatTransportBHELocalAssemblerBHE-impl.h"
