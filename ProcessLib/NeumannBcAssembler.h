/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_NEUMANN_BC_ASSEMBLER_H_
#define PROCESS_LIB_NEUMANN_BC_ASSEMBLER_H_

#include <memory>
#include <vector>

#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "LocalAssemblerUtil.h"

namespace ProcessLib
{

class LocalNeumannBcAsmDataInterface
{
public:
    virtual ~LocalNeumannBcAsmDataInterface() = default;

    virtual void assemble(const double t, std::vector<double>& local_b_data) = 0;
};

template<typename ShapeFunction,
         typename IntegrationMethod,
         unsigned GlobalDim>
class LocalNeumannBcAsmData : public LocalNeumannBcAsmDataInterface
{
public:
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    /// The neumann_bc_value factor is directly integrated into the local
    /// element matrix.
    LocalNeumannBcAsmData(
            MeshLib::Element const& e,
            std::size_t const local_matrix_size,
            unsigned const integration_order,
            std::function<double (MeshLib::Element const&)> const& value_lookup)
        : _shape_matrices(
              initShapeMatrices<ShapeFunction, ShapeMatricesType, IntegrationMethod, GlobalDim>(
                              e, integration_order))
        , _neumann_bc_value(value_lookup(e))
        , _local_matrix_size(local_matrix_size)
        , _integration_order(integration_order)
    {}

    void assemble(const double t, std::vector<double>& local_b_data) override
    {
        (void) t; // TODO time-dependent Neumann BCs

        auto local_b = setupLocalVector(local_b_data, _local_matrix_size);

        IntegrationMethod integration_method(_integration_order);
        std::size_t const n_integration_points = integration_method.getNPoints();

        for (std::size_t ip(0); ip < n_integration_points; ip++) {
            auto const& sm = _shape_matrices[ip];
            auto const& wp = integration_method.getWeightedPoint(ip);
            local_b.noalias() += sm.N * _neumann_bc_value
                        * sm.detJ * wp.getWeight();
        }
    }

private:
    std::vector<ShapeMatrices> _shape_matrices;
    double _neumann_bc_value;

    std::size_t const _local_matrix_size;
    unsigned const _integration_order = 2;
};

}   // namespace ProcessLib

#endif  // PROCESS_LIB_NEUMANN_BC_ASSEMBLER_H_
