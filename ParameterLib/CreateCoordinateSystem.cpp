/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on January 25, 2024, 10:49 AM
 */

#include "CreateCoordinateSystem.h"

#include "BaseLib/ConfigTree.h"
#include "Parameter.h"
#include "Utils.h"

namespace ParameterLib
{
struct ParameterBase;
struct CoordinateSystem;

std::optional<ParameterLib::CoordinateSystem> createCoordinateSystem(
    std::optional<BaseLib::ConfigTree> const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    if (!config)
    {
        return {};
    }

    DBUG("Reading coordinate system configuration.");

    //
    // Fetch the first basis vector; its length defines the dimension.
    //
    auto const& basis_vector_0 = ParameterLib::findParameter<double>(
        *config,
        //! \ogs_file_param_special{prj__local_coordinate_system__basis_vector_0}
        "basis_vector_0", parameters, 0 /* any dimension */);
    int const dimension = basis_vector_0.getNumberOfGlobalComponents();

    // check dimension
    if (dimension != 2 && dimension != 3)
    {
        OGS_FATAL(
            "Basis vector parameter '{:s}' must have two or three components, "
            "but it has {:d}.",
            basis_vector_0.name, dimension);
    }

    //
    // Fetch the second basis vector, which must be of the same dimension as the
    // first one.
    //
    auto const& basis_vector_1 = ParameterLib::findParameter<double>(
        *config,
        //! \ogs_file_param_special{prj__local_coordinate_system__basis_vector_1}
        "basis_vector_1", parameters, dimension);

    //
    // For two dimensions, we are done; construct coordinate system;
    //
    if (dimension == 2)
    {
        return ParameterLib::CoordinateSystem{basis_vector_0, basis_vector_1};
    }

    //
    // Parse the third vector, for three dimensions.
    //
    auto const& basis_vector_2 = ParameterLib::findParameter<double>(
        *config,
        //! \ogs_file_param_special{prj__local_coordinate_system__basis_vector_2}
        "basis_vector_2", parameters, dimension);
    return ParameterLib::CoordinateSystem{basis_vector_0, basis_vector_1,
                                          basis_vector_2};
}
}  // namespace ParameterLib
