/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/ComponentTransport/ComponentTransportFEM.h"
#include "ProcessLib/VectorMatrixAssembler.h"

namespace ProcessLib::ComponentTransport
{

struct AssembledMatrixCache final
{
    AssembledMatrixCache(bool const is_linear,
                         bool const use_monolithic_scheme);

    void assemble(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& x_prev, int const process_id,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
        std::vector<
            std::reference_wrapper<NumLib::LocalToGlobalIndexMap>> const&
            dof_tables,
        VectorMatrixAssembler& global_assembler,
        std::vector<
            std::unique_ptr<ComponentTransportLocalAssemblerInterface>> const&
            local_assemblers,
        std::vector<std::size_t> const& active_element_ids);

    bool isLinear() const { return is_linear_; }

private:
    bool const is_linear_;

    std::unique_ptr<GlobalMatrix> M_{};
    std::unique_ptr<GlobalMatrix> K_{};
    std::unique_ptr<GlobalVector> b_{};
};
}  // namespace ProcessLib::ComponentTransport
