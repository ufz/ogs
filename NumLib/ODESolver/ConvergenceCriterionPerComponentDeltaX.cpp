/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConvergenceCriterionPerComponentDeltaX.h"

#include <limits>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Logging.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"

namespace NumLib
{
ConvergenceCriterionPerComponentDeltaX::ConvergenceCriterionPerComponentDeltaX(
    std::vector<double>&& absolute_tolerances,
    std::vector<double>&& relative_tolerances,
    const MathLib::VecNormType norm_type)
    : ConvergenceCriterionPerComponent(norm_type),
      _abstols(std::move(absolute_tolerances)),
      _reltols(std::move(relative_tolerances))
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

void ConvergenceCriterionPerComponentDeltaX::checkDeltaX(
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

        bool const satisfied_abs = error_dx < _abstols[global_component];
        bool const satisfied_rel = checkRelativeTolerance(
            _reltols[global_component], error_dx, norm_x);

        _satisfied = _satisfied && (satisfied_abs || satisfied_rel);
    }
}

void ConvergenceCriterionPerComponentDeltaX::setDOFTable(
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

    auto const norm_type = MathLib::convertStringToVecNormType(norm_type_str);

    if (norm_type == MathLib::VecNormType::INVALID)
    {
        OGS_FATAL("Unknown vector norm type `{:s}'.", norm_type_str);
    }

    return std::make_unique<ConvergenceCriterionPerComponentDeltaX>(
        std::move(*abstols), std::move(*reltols), norm_type);
}

}  // namespace NumLib
