/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConvergenceCriterionPerComponentResidual.h"
#include "BaseLib/Logging.h"

#include "BaseLib/ConfigTree.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"

namespace NumLib
{
ConvergenceCriterionPerComponentResidual::
    ConvergenceCriterionPerComponentResidual(
        std::vector<double>&& absolute_tolerances,
        std::vector<double>&& relative_tolerances,
        const MathLib::VecNormType norm_type)
    : ConvergenceCriterionPerComponent(norm_type),
      abstols_(std::move(absolute_tolerances)),
      reltols_(std::move(relative_tolerances)),
      residual_norms_0_(abstols_.size())
{
    if (abstols_.size() != reltols_.size())
    {
        OGS_FATAL(
            "The number of absolute and relative tolerances given must be the "
            "same.");
    }

    if (abstols_.empty())
    {
        OGS_FATAL("The given tolerances vector is empty.");
    }
}


void ConvergenceCriterionPerComponentResidual::checkDeltaX(
    const GlobalVector& minus_delta_x, GlobalVector const& x)
{
    if ((!dof_table_) || (!mesh_))
    {
        OGS_FATAL("D.o.f. table or mesh have not been set.");
    }

    for (unsigned global_component = 0; global_component < abstols_.size();
         ++global_component)
    {
        // TODO short cut if tol <= 0.0
        auto error_dx = norm(minus_delta_x, global_component, norm_type_,
                             *dof_table_, *mesh_);
        auto norm_x =
            norm(x, global_component, norm_type_, *dof_table_, *mesh_);

        INFO(
            "Convergence criterion, component {:d}: |dx|={:.4e}, |x|={:.4e}, "
            "|dx|/|x|={:.4e}",
            global_component, error_dx, norm_x,
            (norm_x == 0. ? std::numeric_limits<double>::quiet_NaN()
                          : (error_dx / norm_x)));
    }
}


void ConvergenceCriterionPerComponentResidual::checkResidual(
    const GlobalVector& residual)
{
    if ((!dof_table_) || (!mesh_))
    {
        OGS_FATAL("D.o.f. table or mesh have not been set.");
    }

    bool satisfied_abs = true;
    // Make sure that in the first iteration the relative residual tolerance is
    // not satisfied.
    bool satisfied_rel = !is_first_iteration_;

    for (unsigned global_component = 0; global_component < abstols_.size();
         ++global_component)
    {
        // TODO short cut if tol <= 0.0
        auto norm_res = norm(residual, global_component, norm_type_,
                             *dof_table_, *mesh_);

        if (is_first_iteration_) {
            INFO("Convergence criterion, component {:d}: |r0|={:.4e}",
                 global_component, norm_res);
            residual_norms_0_[global_component] = norm_res;
        } else {
            auto const norm_res0 = residual_norms_0_[global_component];
            INFO(
                "Convergence criterion, component {:d}: |r|={:.4e}, "
                "|r0|={:.4e}, "
                "|r|/|r0|={:.4e}",
                global_component, norm_res, norm_res0,
                (norm_res0 == 0. ? std::numeric_limits<double>::quiet_NaN()
                                 : (norm_res / norm_res0)));
        }

        satisfied_abs = satisfied_abs && norm_res < abstols_[global_component];
        satisfied_rel =
            satisfied_rel &&
            checkRelativeTolerance(reltols_[global_component], norm_res,
                                   residual_norms_0_[global_component]);
    }

    satisfied_ = satisfied_ && (satisfied_abs || satisfied_rel);
}

void ConvergenceCriterionPerComponentResidual::setDOFTable(
    const LocalToGlobalIndexMap& dof_table, MeshLib::Mesh const& mesh)
{
    dof_table_ = &dof_table;
    mesh_ = &mesh;

    if (dof_table_->getNumberOfComponents() !=
        static_cast<int>(abstols_.size()))
    {
        OGS_FATAL(
            "The number of components in the DOF table and the number of "
            "tolerances given do not match.");
    }
}

std::unique_ptr<ConvergenceCriterionPerComponentResidual>
createConvergenceCriterionPerComponentResidual(
    const BaseLib::ConfigTree& config)
{
    //! \ogs_file_param{prj__time_loop__processes__process__convergence_criterion__type}
    config.checkConfigParameter("type", "PerComponentResidual");

    auto abstols =
        //! \ogs_file_param{prj__time_loop__processes__process__convergence_criterion__PerComponentResidual__abstols}
        config.getConfigParameterOptional<std::vector<double>>("abstols");
    auto reltols =
        //! \ogs_file_param{prj__time_loop__processes__process__convergence_criterion__PerComponentResidual__reltols}
        config.getConfigParameterOptional<std::vector<double>>("reltols");
    auto const norm_type_str =
        //! \ogs_file_param{prj__time_loop__processes__process__convergence_criterion__PerComponentResidual__norm_type}
        config.getConfigParameter<std::string>("norm_type");

    if ((!abstols) && (!reltols))
    {
        OGS_FATAL(
            "At least one of absolute or relative tolerance has to be "
            "specified.");
    }
    if (!abstols) {
        abstols = std::vector<double>(reltols->size());
    } else if (!reltols) {
        reltols = std::vector<double>(abstols->size());
    }

    auto const norm_type = MathLib::convertStringToVecNormType(norm_type_str);

    if (norm_type == MathLib::VecNormType::INVALID)
    {
        OGS_FATAL("Unknown vector norm type `{:s}'.", norm_type_str);
    }

    return std::make_unique<ConvergenceCriterionPerComponentResidual>(
        std::move(*abstols), std::move(*reltols), norm_type);
}

}  // namespace NumLib
