// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "NewtonStepStrategy.h"

#include "MathLib/LinAlg/LinAlg.h"
#include "NonlinearSystem.h"
#include "NumLib/DOF/GlobalMatrixProviders.h"

namespace NumLib
{
std::optional<double> NewtonStepContext::evaluateMerit(
    GlobalVector& x_trial) const
{
    // Copy pointers, replace the one for the given process id.
    std::vector<GlobalVector*> x_trials{x};
    x_trials[process_id] = &x_trial;

    if (!assembledOnAllRanks([&]
                             { sys.assemble(x_trials, x_prev, process_id); },
                             "Residual evaluation at a trial point failed.",
                             spdlog::level::warn))
    {
        return std::nullopt;
    }

    auto& res = NumLib::GlobalVectorProvider::provider.getVector(x_trial);
    auto& dirichlet_free_res =
        NumLib::GlobalVectorProvider::provider.getVector(x_trial);
    sys.getResidual(x_trial, *x_prev[process_id], res);
    if (r_neq != nullptr)
    {
        MathLib::LinAlg::axpy(res, -1, *r_neq);
    }
    // x_trial is locally accessible, since the assembly made it so.
    sys.getDirichletFreeResidual(x_trial, res, dirichlet_free_res);
    double const merit =
        MathLib::LinAlg::norm(dirichlet_free_res, MathLib::VecNormType::NORM2);
    NumLib::GlobalVectorProvider::provider.releaseVector(dirichlet_free_res);
    NumLib::GlobalVectorProvider::provider.releaseVector(res);
    return merit;
}
}  // namespace NumLib
