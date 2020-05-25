/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TimeDependentHeterogeneousParameter.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "Utils.h"

namespace ParameterLib
{
TimeDependentHeterogeneousParameter::TimeDependentHeterogeneousParameter(
    std::string name,
    std::vector<PairTimeParameterName>
        time_parameter_name_mapping)
    : Parameter<double>(std::move(name), nullptr),
      time_parameter_name_mapping_(time_parameter_name_mapping)
{
}

int TimeDependentHeterogeneousParameter::getNumberOfComponents() const
{
    return time_parameter_mapping_[0].second->getNumberOfComponents();
}

bool TimeDependentHeterogeneousParameter::isTimeDependent() const
{
    return true;
}

std::vector<double> TimeDependentHeterogeneousParameter::operator()(
    double const t, SpatialPosition const& pos) const
{
    // No local coordinate transformation here, which might happen twice
    // otherwise.
    assert(!this->coordinate_system_ ||
           "Coordinate system not expected to be set for curve scaled "
           "parameters.");
    if (t < time_parameter_mapping_[0].first)
    {
        return time_parameter_mapping_[0].second->operator()(t, pos);
    }
    if (time_parameter_mapping_.back().first <= t)
    {
        return time_parameter_mapping_.back().second->operator()(t, pos);
    }
    std::size_t k(1);
    for (; k < time_parameter_mapping_.size(); ++k)
    {
        if (time_parameter_mapping_[k - 1].first <= t &&
            t < time_parameter_mapping_[k].first)
        {
            break;
        }
    }
    auto const t0 = time_parameter_mapping_[k - 1].first;
    auto const t1 = time_parameter_mapping_[k].first;
    auto const alpha = (t - t0) / (t1 - t0);

    auto r0 = time_parameter_mapping_[k - 1].second->operator()(t, pos);
    std::transform(r0.begin(), r0.end(), r0.begin(),
                   [alpha](auto& v) { return (1 - alpha) * v; });
    auto r1 = time_parameter_mapping_[k].second->operator()(t, pos);
    std::transform(r1.begin(), r1.end(), r1.begin(),
                   [alpha](auto& v) { return alpha * v; });
    std::transform(r0.begin(), r0.end(), r1.begin(), r0.begin(),
                   [](auto& v0, auto const& v1) { return v0 + v1; });
    return r0;
}

void TimeDependentHeterogeneousParameter::initialize(
    std::vector<std::unique_ptr<ParameterBase>> const& parameters)
{
    DBUG("TimeDependentHeterogeneousParameter init {:d} time series entries.",
         time_parameter_name_mapping_.size());
    for (auto const& time_parameter_map : time_parameter_name_mapping_)
    {
        auto parameter =
            &findParameter<double>(time_parameter_map.second, parameters, 0);
        time_parameter_mapping_.emplace_back(time_parameter_map.first,
                                             parameter);
    }

    // check that all parameters have the same number of components
    auto const n = time_parameter_mapping_[0].second->getNumberOfComponents();
    if (!std::all_of(time_parameter_mapping_.begin(),
                     time_parameter_mapping_.end(),
                     [n](auto const p) {
                         return n == p.second->getNumberOfComponents();
                     }))
    {
        OGS_FATAL(
            "All referenced parameters in time dependent heterogeneous "
            "parameter '{:s}' have to have the same number of components.",
            name);
    }
}

std::unique_ptr<ParameterBase> createTimeDependentHeterogeneousParameter(
    std::string const& name, BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__parameters__parameter__type}
    config.checkConfigParameter("type", "TimeDependentHeterogeneousParameter");
    auto const time_series_config =
        //! \ogs_file_param{prj__parameters__parameter__TimeDependentHeterogeneousParameter__time_series}
        config.getConfigSubtree("time_series");

    std::vector<TimeDependentHeterogeneousParameter::PairTimeParameterName>
        time_series;
    //! \ogs_file_param{prj__parameters__parameter__TimeDependentHeterogeneousParameter__time_series__pair}
    for (auto const p : time_series_config.getConfigSubtreeList("pair"))
    {
        //! \ogs_file_param{prj__parameters__parameter__TimeDependentHeterogeneousParameter__time_series__pair__time}
        auto time = p.getConfigParameter<double>("time");
        auto parameter_name =
        //! \ogs_file_param{prj__parameters__parameter__TimeDependentHeterogeneousParameter__time_series__pair__parameter_name}
            p.getConfigParameter<std::string>("parameter_name");
        time_series.emplace_back(time, parameter_name);
    }

    if (time_series.empty())
    {
        OGS_FATAL(
            "Time dependent heterogeneous parameter '{:s}' doesn't contain "
            "necessary time series data.",
            name);
    }

    if (!std::is_sorted(
            time_series.begin(), time_series.end(),
            [](TimeDependentHeterogeneousParameter::PairTimeParameterName const&
                   p0,
               TimeDependentHeterogeneousParameter::PairTimeParameterName const&
                   p1) { return p0.first < p1.first; }))
    {
        OGS_FATAL(
            "The points in time in the time series '{:s}' aren't in ascending "
            "order.",
            name);
    }

    return std::make_unique<TimeDependentHeterogeneousParameter>(
        name, std::move(time_series));
}

}  // namespace ParameterLib
