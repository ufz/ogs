/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "LevelSetFunction.h"

#include <boost/math/special_functions/sign.hpp>
#include <numeric>

#include "BaseLib/Algorithm.h"
#include "BranchProperty.h"
#include "FractureProperty.h"
#include "JunctionProperty.h"

namespace
{
// Heaviside step function
inline double Heaviside(bool const v)
{
    return v ? 0.5 : -0.5;
}

}  // namespace

namespace ProcessLib
{
namespace LIE
{
bool levelsetFracture(FractureProperty const& frac, Eigen::Vector3d const& x)
{
    return frac.normal_vector.dot(x - frac.point_on_fracture) > 0;
}

bool levelsetBranch(BranchProperty const& branch, Eigen::Vector3d const& x)
{
    return branch.normal_vector_branch.dot(x - branch.coords) > 0;
}

std::vector<double> uGlobalEnrichments(
    std::vector<FractureProperty*> const& frac_props,
    std::vector<JunctionProperty*> const& junction_props,
    std::unordered_map<int, int> const& fracID_to_local,
    Eigen::Vector3d const& x)
{
    // pre-calculate levelsets for all fractures
    std::vector<bool> levelsets(frac_props.size());
    for (std::size_t i = 0; i < frac_props.size(); i++)
    {
        levelsets[i] = levelsetFracture(*frac_props[i], x);
    }

    std::vector<double> enrichments(frac_props.size() + junction_props.size());
    // fractures possibly with branches
    for (std::size_t i = 0; i < frac_props.size(); i++)
    {
        auto const* frac = frac_props[i];
        enrichments[i] = Heaviside(std::accumulate(
            cbegin(frac->branches_slave), cend(frac->branches_slave),
            levelsets[i], [&](bool const enrich, auto const& branch) {
                return enrich & levelsetBranch(branch, x);
            }));
    }

    // junctions
    for (std::size_t i = 0; i < junction_props.size(); i++)
    {
        auto const* junction = junction_props[i];
        auto fid1 = fracID_to_local.at(junction->fracture_ids[0]);
        auto fid2 = fracID_to_local.at(junction->fracture_ids[1]);
        bool const enrich = levelsets[fid1] & levelsets[fid2];
        enrichments[i + frac_props.size()] = Heaviside(enrich);
    }

    return enrichments;
}

std::vector<double> duGlobalEnrichments(
    std::size_t this_frac_id,
    std::vector<FractureProperty*> const& frac_props,
    std::vector<JunctionProperty*> const& junction_props,
    std::unordered_map<int, int> const& fracID_to_local,
    Eigen::Vector3d const& x)
{
    auto this_frac_local_index = fracID_to_local.at(this_frac_id);
    auto const& this_frac = *frac_props[this_frac_local_index];
    // pre-calculate levelsets for all fractures
    std::vector<bool> levelsets(frac_props.size());
    for (std::size_t i = 0; i < frac_props.size(); i++)
    {
        levelsets[i] = levelsetFracture(*frac_props[i], x);
    }

    std::vector<double> enrichments(frac_props.size() + junction_props.size());
    enrichments[this_frac_local_index] = 1.0;

    // fractures possibly with branches
    if (frac_props.size() > 1)
    {
        for (auto const& branch : this_frac.branches_master)
        {
            if (branch.master_fracture_id != this_frac.fracture_id)
            {
                continue;
            }

            if (fracID_to_local.find(branch.slave_fracture_id) ==
                fracID_to_local.end())
            {
                continue;
            }

            double sign = boost::math::sign(
                this_frac.normal_vector.dot(branch.normal_vector_branch));
            auto slave_fid = fracID_to_local.at(branch.slave_fracture_id);
            double const enrich = levelsets[slave_fid] ? 1. : 0.;
            enrichments[slave_fid] = sign * enrich;
        }
    }

    // junctions
    for (unsigned i = 0; i < junction_props.size(); i++)
    {
        auto const* junction = junction_props[i];
        if (!BaseLib::contains(junction->fracture_ids, this_frac.fracture_id))
        {
            continue;
        }

        auto another_frac_id =
            (junction->fracture_ids[0] == this_frac.fracture_id)
                ? junction->fracture_ids[1]
                : junction->fracture_ids[0];
        auto fid = fracID_to_local.at(another_frac_id);
        double const enrich = levelsets[fid] ? 1. : 0.;
        enrichments[i + frac_props.size()] = enrich;
    }

    return enrichments;
}

}  // namespace LIE
}  // namespace ProcessLib
