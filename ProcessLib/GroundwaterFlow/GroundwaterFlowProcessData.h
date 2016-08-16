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

template <typename T>
struct Parameter;

namespace GroundwaterFlow
{

struct GroundwaterFlowProcessData
{
    GroundwaterFlowProcessData(
        Parameter<double> const& hydraulic_conductivity_)
        : hydraulic_conductivity(hydraulic_conductivity_)
    {}

    GroundwaterFlowProcessData(GroundwaterFlowProcessData&& other)
        : hydraulic_conductivity(other.hydraulic_conductivity)
    {}

    //! Copies are forbidden.
    GroundwaterFlowProcessData(GroundwaterFlowProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(GroundwaterFlowProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(GroundwaterFlowProcessData&&) = delete;

    Parameter<double> const& hydraulic_conductivity;
};

} // namespace GroundwaterFlow
} // namespace ProcessLib

#endif // PROCESSLIB_GROUNDWATERFLOW_GROUNDWATERFLOWPROCESSDATA_H
