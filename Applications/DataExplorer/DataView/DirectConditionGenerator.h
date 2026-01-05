// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
