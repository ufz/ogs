// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "VtkCompositeFilter.h"
#include "GeoLib/Point.h"

#include <vector>

/// @brief This filter displays the points/nodes given in the index field as spheres.
class VtkCompositeNodeSelectionFilter : public VtkCompositeFilter
{
public:
    explicit VtkCompositeNodeSelectionFilter(vtkAlgorithm* inputAlgorithm);
    ~VtkCompositeNodeSelectionFilter() override;

    void init() override;

    /// Sets the point indices to be highlighted
    void setSelectionArray(const std::vector<unsigned> &point_indeces);

private:
    std::vector<GeoLib::Point*> _selection;
};
