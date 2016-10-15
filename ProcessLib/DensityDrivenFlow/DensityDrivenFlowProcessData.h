/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_DensityDrivenFlow_DensityDrivenFlowPROCESSDATA_H
#define PROCESSLIB_DensityDrivenFlow_DensityDrivenFlowPROCESSDATA_H

namespace MeshLib
{
    class Element;
}


namespace ProcessLib
{

template <typename ReturnType>
struct Parameter;

namespace DensityDrivenFlow
{

// copy that and change names first still only thermal conductivity, later on more parameters will be added.
struct DensityDrivenFlowProcessData
{
    DensityDrivenFlowProcessData(
            ProcessLib::Parameter<double> const&
            thermal_conductivity_
            )
        : thermal_conductivity(thermal_conductivity_)
    {}

    DensityDrivenFlowProcessData(DensityDrivenFlowProcessData&& other)
        : thermal_conductivity(other.thermal_conductivity)
    {}

    //! Copies are forbidden.
    DensityDrivenFlowProcessData(DensityDrivenFlowProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(DensityDrivenFlowProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(DensityDrivenFlowProcessData&&) = delete;

    Parameter<double> const& thermal_conductivity;
};

} // namespace DensityDrivenFlow
} // namespace ProcessLib

#endif // PROCESSLIB_DensityDrivenFlow_DensityDrivenFlowPROCESSDATA_H
