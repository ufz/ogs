/**
 * \file
 * \author Karsten Rink
 * \date   2013-04-16
 * \brief  Definition of the VtkCompositeNodeSelectionFilter class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "VtkCompositeFilter.h"
#include "Point.h"

#include <vector>

/// @brief This filter displays the points/nodes given in the index field as spheres.
class VtkCompositeNodeSelectionFilter : public VtkCompositeFilter
{
public:
    VtkCompositeNodeSelectionFilter(vtkAlgorithm* inputAlgorithm);
    ~VtkCompositeNodeSelectionFilter() override;

    void init() override;

    /// Sets the point indices to be highlighted
    void setSelectionArray(const std::vector<unsigned> &point_indeces);

private:
    std::vector<GeoLib::Point*> _selection;
};
