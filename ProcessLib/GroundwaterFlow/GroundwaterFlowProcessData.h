/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_GROUNDWATERFLOW_GROUNDWATERFLOWPROCESSDATA_H
#define PROCESSLIB_GROUNDWATERFLOW_GROUNDWATERFLOWPROCESSDATA_H

namespace MeshLib
{
    class Element;
}


namespace ProcessLib
{

template <typename ReturnType, typename... Args>
struct Parameter;

namespace GroundwaterFlow
{

struct GroundwaterFlowProcessData
{
    GroundwaterFlowProcessData(
            ProcessLib::Parameter<double, MeshLib::Element const&> const&
            hydraulic_conductivity_
            )
        : hydraulic_conductivity(hydraulic_conductivity_)
    {}

    GroundwaterFlowProcessData(GroundwaterFlowProcessData&& other)
        : hydraulic_conductivity(other.hydraulic_conductivity)
    {}

    //! Copies are forbidden.
    GroundwaterFlowProcessData(GroundwaterFlowProcessData const&) = delete;

    Parameter<double, MeshLib::Element const&> const& hydraulic_conductivity;
};

} // namespace GroundwaterFlow
} // namespace ProcessLib

#endif // PROCESSLIB_GROUNDWATERFLOW_GROUNDWATERFLOWPROCESSDATA_H
