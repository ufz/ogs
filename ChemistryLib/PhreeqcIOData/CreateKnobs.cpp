/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <boost/optional/optional.hpp>

#include "BaseLib/ConfigTree.h"
#include "Knobs.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
std::unique_ptr<Knobs> createKnobs(
    boost::optional<BaseLib::ConfigTree> const& config)
{
    if (!config)
    {
        return nullptr;
    }

    auto const max_iterations =
        //! \ogs_file_param{prj__chemical_system__knobs__max_iterations}
        config->getConfigParameter<int>("max_iter", 100);

    auto const relative_convergence_tolerance =
        //! \ogs_file_param{prj__chemical_system__knobs__relative_convergence_tolerance}
        config->getConfigParameter<double>("relative_convergence_tolerance",
                                           1e-12);

    auto const tolerance =
        //! \ogs_file_param{prj__chemical_system__knobs__tolerance}
        config->getConfigParameter<double>("tolerance", 1e-15);

    auto const step_size =
        //! \ogs_file_param{prj__chemical_system__knobs__tolerance}
        config->getConfigParameter<int>("step_size", 100);

    auto const scaling =
        //! \ogs_file_param{prj__chemical_system__knobs__scaling}
        config->getConfigParameter<bool>("scaling", false);

    return std::make_unique<Knobs>(max_iterations,
                                   relative_convergence_tolerance, tolerance,
                                   step_size, scaling);
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
