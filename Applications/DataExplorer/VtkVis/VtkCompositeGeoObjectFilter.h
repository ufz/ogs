/**
 * \file
 * \author Karsten Rink
 * \date   2011-12-02
 * \brief  Definition of the VtkCompositeGeoObjectFilter class.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "VtkCompositeFilter.h"
#include "GeoLib/GeoType.h"

class vtkThreshold;

/// @brief Highlights a single GeoObject
class VtkCompositeGeoObjectFilter : public VtkCompositeFilter
{
public:
    explicit VtkCompositeGeoObjectFilter(vtkAlgorithm* inputAlgorithm);
    ~VtkCompositeGeoObjectFilter() override;

    void init() override;

    /// @brief Sets user properties.
    void SetUserProperty(QString name, QVariant value) override
    {
        Q_UNUSED(name);
        Q_UNUSED(value);
    }

    void SetIndex(std::size_t idx);

private:
    GeoLib::GEOTYPE _type;
    vtkThreshold* _threshold;
};
