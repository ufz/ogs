/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "PythonSourceTerm.h"
#include "Utils/BcAndStLocalAssemblerImpl.h"

namespace ProcessLib
{
namespace SourceTerms
{
namespace Python
{
template <typename ShapeFunction, typename LowerOrderShapeFunction,
          typename IntegrationMethod, int GlobalDim>
class PythonSourceTermLocalAssembler final
    : public PythonSourceTermLocalAssemblerInterface
{
    using LocAsmImpl = ProcessLib::BoundaryConditionAndSourceTerm::Python::
        BcAndStLocalAssemblerImpl<PythonStData, ShapeFunction,
                                  LowerOrderShapeFunction, IntegrationMethod,
                                  GlobalDim>;
    using Traits = typename LocAsmImpl::Traits;

public:
    PythonSourceTermLocalAssembler(MeshLib::Element const& e,
                                   std::size_t const /*local_matrix_size*/,
                                   bool const is_axially_symmetric,
                                   unsigned const integration_order,
                                   PythonStData const& data)
        : impl_{e, is_axially_symmetric, integration_order, data}
    {
    }

    void assemble(std::size_t const source_term_element_id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_source_term,
                  double const t, const GlobalVector& x, GlobalVector& b,
                  GlobalMatrix* const Jac) override
    {
        impl_.assemble(source_term_element_id, dof_table_source_term, t, x, b,
                       Jac);
    }

private:
    LocAsmImpl const impl_;
};

}  // namespace Python
}  // namespace SourceTerms
}  // namespace ProcessLib
