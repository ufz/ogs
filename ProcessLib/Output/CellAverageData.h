// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MeshLib/Mesh.h"
#include "MeshLib/PropertyVector.h"

namespace ProcessLib
{
struct CellAverageData
{
    explicit CellAverageData(MeshLib::Mesh& mesh) : mesh_{mesh} {}

    MeshLib::PropertyVector<double>& getOrCreatePropertyVector(
        std::string const& name, unsigned const num_comp);

private:
    MeshLib::Mesh const& mesh_;
    std::map<std::string, MeshLib::PropertyVector<double>*> cell_averages_;
};
}  // namespace ProcessLib
