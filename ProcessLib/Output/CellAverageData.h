/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
