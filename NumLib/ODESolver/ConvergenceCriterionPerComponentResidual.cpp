/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConvergenceCriterionPerComponentResidual.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Logging.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"

namespace NumLib
{
ConvergenceCriterionPerComponentResidual::
    ConvergenceCriterionPerComponentResidual(
        std::vector<double>&& absolute_tolerances,
        std::vector<double>&& relative_tolerances,
        std::vector<double>&& damping_alpha,
        bool damping_alpha_switch,
        const MathLib::VecNormType norm_type)
    : ConvergenceCriterionPerComponent(norm_type),
      _abstols(std::move(absolute_tolerances)),
      _reltols(std::move(relative_tolerances)),
      _residual_norms_0(_abstols.size()),
      _damping_alpha(std::move(damping_alpha)),
      _damping_alpha_switch(damping_alpha_switch)
{
    if (_abstols.size() != _reltols.size())
    {
        OGS_FATAL(
            "The number of absolute and relative tolerances given must be the "
            "same.");
    }

    if (_abstols.empty())
    {
        OGS_FATAL("The given tolerances vector is empty.");
    }
}

void ConvergenceCriterionPerComponentResidual::checkDeltaX(
    const GlobalVector& minus_delta_x, GlobalVector const& x)
{
    if ((!_dof_table) || (!_mesh))
    {
        OGS_FATAL("D.o.f. table or mesh have not been set.");
    }

    for (unsigned global_component = 0; global_component < _abstols.size();
         ++global_component)
    {
        // TODO short cut if tol <= 0.0
        auto error_dx = norm(minus_delta_x, global_component, _norm_type,
                             *_dof_table, *_mesh);
        auto norm_x =
            norm(x, global_component, _norm_type, *_dof_table, *_mesh);

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
    if ((!_dof_table) || (!_mesh))
    {
        OGS_FATAL("D.o.f. table or mesh have not been set.");
    }

    for (unsigned global_component = 0; global_component < _abstols.size();
         ++global_component)
    {
        // TODO short cut if tol <= 0.0
        auto norm_res =
            norm(residual, global_component, _norm_type, *_dof_table, *_mesh);

        if (_is_first_iteration)
        {
            INFO("Convergence criterion, component {:d}: |r0|={:.4e}",
                 global_component, norm_res);
            _residual_norms_0[global_component] = norm_res;
        }
        else
        {
            auto const norm_res0 = _residual_norms_0[global_component];
            INFO(
                "Convergence criterion, component {:d}: |r|={:.4e}, "
                "|r0|={:.4e}, |r|/|r0|={:.4e}",
                global_component, norm_res, norm_res0,
                (norm_res0 == 0. ? std::numeric_limits<double>::quiet_NaN()
                                 : (norm_res / norm_res0)));
        }

        bool const satisfied_abs = norm_res < _abstols[global_component];
        // Make sure that in the first iteration the relative residual tolerance
        // is not satisfied.
        bool const satisfied_rel =
            !_is_first_iteration &&
            checkRelativeTolerance(_reltols[global_component], norm_res,
                                   _residual_norms_0[global_component]);
        _satisfied = _satisfied && (satisfied_abs || satisfied_rel);
    }
}

double ConvergenceCriterionPerComponentResidual::getDampingFactor(
    const GlobalVector& minus_delta_x, GlobalVector const& x,
    double damping_orig)
{
    if ((!_dof_table) || (!_mesh))
    {
        OGS_FATAL("D.o.f. table or mesh have not been set.");
    }

    double damping_final = 1;
    for (unsigned global_component = 0;
         global_component < _damping_alpha.size();
         ++global_component)
    {
        auto const& ms = _dof_table->getMeshSubset(global_component);
        assert(ms.getMeshID() == _mesh->getID());
        DBUG("Non-negative damping for component: {:d} alpha: {:g}",
             global_component, _damping_alpha[global_component]);
        for (MeshLib::Node const* node : ms.getNodes())
        {
            MeshLib::Location const l{
                _mesh->getID(), MeshLib::MeshItemType::Node, node->getID()};
            auto index = _dof_table->getGlobalIndex(l, global_component);
            auto const damping_tmp =
                damping_orig / std::max(1.0, (minus_delta_x.get(index) *
                                              _damping_alpha[global_component] /
                                              x.get(index)));
            damping_final =
                (damping_tmp < damping_final) ? damping_tmp : damping_final;
        }
    }
    DBUG("Final damping value due to non-negative damping: {:g}",
         damping_final);
    return damping_final;
}

void ConvergenceCriterionPerComponentResidual::setDOFTable(
    const LocalToGlobalIndexMap& dof_table, MeshLib::Mesh const& mesh)
{
    _dof_table = &dof_table;
    _mesh = &mesh;

    if (_dof_table->getNumberOfGlobalComponents() !=
        static_cast<int>(_abstols.size()))
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
    auto damping_alpha =
        //! \ogs_file_param{prj__time_loop__processes__process__convergence_criterion__PerComponentResidual__damping_alpha}
        config.getConfigParameterOptional<std::vector<double>>("damping_alpha");
    auto const norm_type_str =
        //! \ogs_file_param{prj__time_loop__processes__process__convergence_criterion__PerComponentResidual__norm_type}
        config.getConfigParameter<std::string>("norm_type");

    bool damping_alpha_switch = true;
    if ((!abstols) && (!reltols))
    {
        OGS_FATAL(
            "At least one of absolute or relative tolerance has to be "
            "specified.");
    }
    if (!abstols)
    {
        abstols = std::vector<double>(reltols->size());
    }
    else if (!reltols)
    {
        reltols = std::vector<double>(abstols->size());
    }
    if (!damping_alpha)
    {
        damping_alpha = std::vector<double>(abstols->size(), 0.0);
        damping_alpha_switch = false;
    }

    auto const norm_type = MathLib::convertStringToVecNormType(norm_type_str);

    if (norm_type == MathLib::VecNormType::INVALID)
    {
        OGS_FATAL("Unknown vector norm type `{:s}'.", norm_type_str);
    }

    return std::make_unique<ConvergenceCriterionPerComponentResidual>(
        std::move(*abstols), std::move(*reltols), std::move(*damping_alpha),
        damping_alpha_switch, norm_type);
}

}  // namespace NumLib
