/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConvergenceCriterionPerComponentDeltaX.h"
#include <logog/include/logog.hpp>

#include "BaseLib/ConfigTree.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/DOF/DOFTableUtil.h"

namespace NumLib
{
ConvergenceCriterionPerComponentDeltaX::ConvergenceCriterionPerComponentDeltaX(
    std::vector<double>&& absolute_tolerances,
    std::vector<double>&& relative_tolerances,
    MathLib::VecNormType norm_type)
    : _abstols(std::move(absolute_tolerances)),
      _reltols(std::move(relative_tolerances)),
      _norm_type(norm_type)
{
    if (_abstols.size() != _reltols.size())
        OGS_FATAL(
            "The number of absolute and relative tolerances given must be the "
            "same.");

    if (_abstols.empty())
        OGS_FATAL("The given tolerances vector is empty.");
}

void ConvergenceCriterionPerComponentDeltaX::checkDeltaX(
    const GlobalVector& minus_delta_x, GlobalVector const& x)
{
    if ((!_dof_table) || (!_mesh))
        OGS_FATAL("D.o.f. table or mesh have not been set.");

    bool satisfied_abs = true;
    bool satisfied_rel = true;

    for (unsigned global_component = 0; global_component < _abstols.size();
         ++global_component)
    {
        // TODO short cut if tol <= 0.0
        auto error_dx = norm(minus_delta_x, global_component, _norm_type,
                             *_dof_table, *_mesh);
        auto norm_x =
            norm(x, global_component, _norm_type, *_dof_table, *_mesh);

        INFO(
            "Convergence criterion, component %u: |dx|=%.4e, |x|=%.4e, "
            "|dx|/|x|=%.4e",
            error_dx, global_component, norm_x, error_dx / norm_x);

        satisfied_abs = satisfied_abs && error_dx < _abstols[global_component];
        satisfied_rel =
            satisfied_rel && checkRelativeTolerance(_reltols[global_component],
                                                    error_dx, norm_x);
    }

    _satisfied = _satisfied && (satisfied_abs || satisfied_rel);
}

void ConvergenceCriterionPerComponentDeltaX::setDOFTable(
    const LocalToGlobalIndexMap& dof_table, MeshLib::Mesh const& mesh)
{
    _dof_table = &dof_table;
    _mesh = &mesh;

    if (_dof_table->getNumberOfComponents() != _abstols.size())
        OGS_FATAL(
            "The number of components in the DOF table and the number of "
            "tolerances given do not match.");
}

std::unique_ptr<ConvergenceCriterionPerComponentDeltaX>
createConvergenceCriterionPerComponentDeltaX(const BaseLib::ConfigTree& config)
{
    //! \ogs_file_param{prj__time_loop__processes__process__convergence_criterion__type}
    config.checkConfigParameter("type", "PerComponentDeltaX");

    auto abstols =
        //! \ogs_file_param{prj__time_loop__processes__process__convergence_criterion__PerComponentDeltaX__abstols}
        config.getConfigParameterOptional<std::vector<double>>("abstols");
    auto reltols =
        //! \ogs_file_param{prj__time_loop__processes__process__convergence_criterion__PerComponentDeltaX__reltols}
        config.getConfigParameterOptional<std::vector<double>>("reltols");
    auto const norm_type_str =
        //! \ogs_file_param{prj__time_loop__processes__process__convergence_criterion__PerComponentDeltaX__norm_type}
        config.getConfigParameter<std::string>("norm_type");

    if ((!abstols) && (!reltols))
        OGS_FATAL(
            "At least one of absolute or relative tolerance has to be "
            "specified.");
    if (!abstols) {
        abstols = std::vector<double>(reltols->size());
    } else if (!reltols) {
        reltols = std::vector<double>(abstols->size());
    }

    auto const norm_type = MathLib::convertStringToVecNormType(norm_type_str);

    if (norm_type == MathLib::VecNormType::INVALID)
        OGS_FATAL("Unknown vector norm type `%s'.", norm_type_str.c_str());

    return std::unique_ptr<ConvergenceCriterionPerComponentDeltaX>(
        new ConvergenceCriterionPerComponentDeltaX(
            std::move(*abstols), std::move(*reltols), norm_type));
}

}  // NumLib
