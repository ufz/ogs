/**
 * \file
 * \author Karsten Rink
 * \date   2012-01-04
 * \brief  Definition of the DirectConditionGenerator class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>
#include <string>

namespace MeshLib {
    class Mesh;
}

class DirectConditionGenerator
{
public:
    DirectConditionGenerator() = default;
    ~DirectConditionGenerator() = default;

    const std::vector<std::pair<std::size_t, double>>& directToSurfaceNodes(
        const MeshLib::Mesh& mesh, const std::string& filename);

    const std::vector<std::pair<std::size_t, double>>&
    directWithSurfaceIntegration(MeshLib::Mesh& mesh,
                                 const std::string& filename, double scaling);

    int writeToFile(const std::string& name) const;

private:
    std::vector< std::pair<std::size_t,double> > _direct_values;

};
