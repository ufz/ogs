/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <iterator>
#include <map>
#include <vector>

#include "BaseLib/Logging.h"
#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"

namespace ProcessLib
{
namespace ComponentTransport
{
struct InterpolationPoint
{
    std::pair<double, double> bounding_points;
    double value;
};

struct Field
{
    Field(std::vector<std::vector<std::size_t>> point_id_groups_,
          std::vector<double> seed_points_, std::string name_, int variable_id_)
        : point_id_groups(std::move(point_id_groups_)),
          seed_points(std::move(seed_points_)),
          name(std::move(name_)),
          variable_id(variable_id_)
    {
    }

    std::pair<double, double> getBoundingSeedPoints(double const value) const;

    std::vector<std::vector<std::size_t>> const point_id_groups;
    std::vector<double> const seed_points;
    std::string const name;
    int const variable_id;
};

struct LookupTable
{
    LookupTable(std::vector<Field> input_fields_,
                std::map<std::string, std::vector<double>>
                    tabular_data_)
        : input_fields(std::move(input_fields_)),
          tabular_data(std::move(tabular_data_))
    {
    }

    void lookup(std::vector<GlobalVector*> const& x,
                std::vector<GlobalVector*> const& x_prev,
                std::size_t const n_nodes) const;

    std::size_t getTableEntryID(std::vector<double> const& entry_input) const;

    std::vector<Field> const input_fields;
    std::map<std::string, std::vector<double>> const tabular_data;
};
}  // namespace ComponentTransport
}  // namespace ProcessLib
