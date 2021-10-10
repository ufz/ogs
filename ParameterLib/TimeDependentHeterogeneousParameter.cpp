/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
      _time_parameter_name_mapping(std::move(time_parameter_name_mapping))
{
}

int TimeDependentHeterogeneousParameter::getNumberOfGlobalComponents() const
{
    return _time_parameter_mapping[0].second->getNumberOfGlobalComponents();
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
    assert(!this->_coordinate_system ||
           "Coordinate system not expected to be set for curve scaled "
           "parameters.");
    if (t < _time_parameter_mapping[0].first)
    {
        return _time_parameter_mapping[0].second->operator()(t, pos);
    }
    if (_time_parameter_mapping.back().first <= t)
    {
        return _time_parameter_mapping.back().second->operator()(t, pos);
    }
    std::size_t k(1);
    for (; k < _time_parameter_mapping.size(); ++k)
    {
        if (_time_parameter_mapping[k - 1].first <= t &&
            t < _time_parameter_mapping[k].first)
        {
            break;
        }
    }
    auto const t0 = _time_parameter_mapping[k - 1].first;
    auto const t1 = _time_parameter_mapping[k].first;
    auto const alpha = (t - t0) / (t1 - t0);

    auto r0 = _time_parameter_mapping[k - 1].second->operator()(t, pos);
    std::transform(r0.begin(), r0.end(), r0.begin(),
                   [alpha](auto const& v) { return (1 - alpha) * v; });
    auto r1 = _time_parameter_mapping[k].second->operator()(t, pos);
    std::transform(r1.begin(), r1.end(), r1.begin(),
                   [alpha](auto const& v) { return alpha * v; });
    std::transform(r0.begin(), r0.end(), r1.begin(), r0.begin(),
                   [](auto const& v0, auto const& v1) { return v0 + v1; });
    return r0;
}

void TimeDependentHeterogeneousParameter::initialize(
    std::vector<std::unique_ptr<ParameterBase>> const& parameters)
{
    DBUG("TimeDependentHeterogeneousParameter init {:d} time series entries.",
         _time_parameter_name_mapping.size());
    for (auto const& time_parameter_map : _time_parameter_name_mapping)
    {
        auto parameter =
            &findParameter<double>(time_parameter_map.second, parameters, 0);
        _time_parameter_mapping.emplace_back(time_parameter_map.first,
                                             parameter);
    }

    // check that all parameters have the same number of components
    auto const n =
        _time_parameter_mapping[0].second->getNumberOfGlobalComponents();
    if (!std::all_of(_time_parameter_mapping.begin(),
                     _time_parameter_mapping.end(),
                     [n](auto const p)
                     { return n == p.second->getNumberOfGlobalComponents(); }))
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
