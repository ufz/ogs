/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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
const unsigned NUM_NODAL_DOF_SOIL = 1;

template <typename ShapeFunction, typename IntegrationMethod, int GlobalDim>
class HeatTransportBHELocalAssemblerSoil
    : public HeatTransportBHELocalAssemblerInterface
{
public:
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    HeatTransportBHELocalAssemblerSoil(
        HeatTransportBHELocalAssemblerSoil const&) = delete;
    HeatTransportBHELocalAssemblerSoil(HeatTransportBHELocalAssemblerSoil&&) =
        delete;

    HeatTransportBHELocalAssemblerSoil(
        MeshLib::Element const& e,
        std::size_t const local_matrix_size,
        std::vector<unsigned> const& dofIndex_to_localIndex,
        bool is_axially_symmetric,
        unsigned const integration_order,
        HeatTransportBHEProcessData& process_data);

    void assemble(double const /*t*/, std::vector<double> const& /*local_x*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_b_data*/) override;

    void assembleWithJacobian(double const /*t*/,
                              std::vector<double> const& /*local_x*/,
                              std::vector<double> const& /*local_xdot*/,
                              const double /*dxdot_dx*/, const double /*dx_dx*/,
                              std::vector<double>& /*local_M_data*/,
                              std::vector<double>& /*local_K_data*/,
                              std::vector<double>& /*local_b_data*/,
                              std::vector<double>& /*local_Jac_data*/) override
    {
        OGS_FATAL(
            "HeatTransportBHELocalAssemblerMatrix: assembly with jacobian is "
            "not "
            "implemented.");
    }

    void preTimestepConcrete(std::vector<double> const& /*local_x*/,
                             double const /*t*/,
                             double const /*delta_t*/) override
    {
    }

    void postTimestepConcrete(std::vector<double> const& /*local_x*/) override;

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
        IntegrationPointDataSoil<ShapeMatricesType>,
        Eigen::aligned_allocator<IntegrationPointDataSoil<ShapeMatricesType>>>
        _ip_data;

    IntegrationMethod const _integration_method;

    std::vector<ShapeMatrices, Eigen::aligned_allocator<ShapeMatrices>>
        _shape_matrices;

    std::size_t const element_id;

    bool const _is_axially_symmetric;

    SecondaryData<typename ShapeMatrices::ShapeType> _secondary_data;
};
}  // namespace HeatTransportBHE
}  // namespace ProcessLib

#include "HeatTransportBHELocalAssemblerSoil_impl.h"
