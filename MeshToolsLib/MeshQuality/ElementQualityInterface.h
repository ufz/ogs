/**
 * \file
 * \author Karsten Rink
 * \date   2015-03-24
 * \brief  Definition of the ElementQualityInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <vector>

#include "BaseLib/Histogram.h"
#include "MeshLib/Mesh.h"
#include "MeshToolsLib/MeshQuality/AngleSkewMetric.h"
#include "MeshToolsLib/MeshQuality/EdgeRatioMetric.h"
#include "MeshToolsLib/MeshQuality/ElementQualityMetric.h"
#include "MeshToolsLib/MeshQuality/ElementSizeMetric.h"
#include "MeshToolsLib/MeshQuality/RadiusEdgeRatioMetric.h"
#include "MeshToolsLib/MeshQuality/SizeDifferenceMetric.h"

namespace MeshToolsLib
{

/**
 * Interface class for handling mesh element quality metrics
 */
class ElementQualityInterface
{
public:
    /// Constructor
    ElementQualityInterface(MeshLib::Mesh const& mesh,
                            MeshLib::MeshQualityType t)
        : _type(t), _mesh(mesh), _quality_tester(nullptr)
    {
        calculateElementQuality(_mesh, _type);
    }

    /// Returns the vector containing a quality measure for each element.
    std::vector<double> const getQualityVector() const
    {
        if (_quality_tester)
            return _quality_tester->getElementQuality();

        std::vector<double> empty_quality_vec(0);
        return empty_quality_vec;
    }

    /// Returns a histogram of the quality vector separated into the given
    /// number of bins. If no number of bins is specified, one will be
    /// calculated based on the Sturges criterium.
    BaseLib::Histogram<double> getHistogram(std::size_t n_bins = 0) const
    {
        if (_quality_tester)
            return _quality_tester->getHistogram(
                static_cast<std::size_t>(n_bins));

        return BaseLib::Histogram<double>{{}};
    }

    /// Writes a histogram of the quality vector to a specified file.
    int writeHistogram(std::string const& file_name,
                       std::size_t n_bins = 0) const
    {
        if (_quality_tester == nullptr)
            return 1;

        BaseLib::Histogram<double> const histogram(
            _quality_tester->getHistogram(n_bins));
        histogram.write(file_name, _mesh.getName(),
                        MeshQualityType2String(_type));
        return 0;
    }

private:
    /// Calculates the quality of each mesh element based on the specified
    /// metric
    void calculateElementQuality(MeshLib::Mesh const& mesh,
                                 MeshLib::MeshQualityType t)
    {
        if (t == MeshLib::MeshQualityType::EDGERATIO)
            _quality_tester =
                std::make_unique<MeshToolsLib::EdgeRatioMetric>(mesh);
        else if (t == MeshLib::MeshQualityType::ELEMENTSIZE)
            _quality_tester =
                std::make_unique<MeshToolsLib::ElementSizeMetric>(mesh);
        else if (t == MeshLib::MeshQualityType::SIZEDIFFERENCE)
            _quality_tester =
                std::make_unique<MeshToolsLib::SizeDifferenceMetric>(mesh);
        else if (t == MeshLib::MeshQualityType::EQUIANGLESKEW)
            _quality_tester =
                std::make_unique<MeshToolsLib::AngleSkewMetric>(mesh);
        else if (t == MeshLib::MeshQualityType::RADIUSEDGERATIO)
            _quality_tester =
                std::make_unique<MeshToolsLib::RadiusEdgeRatioMetric>(mesh);
        else
        {
            ERR("ElementQualityInterface::calculateElementQuality(): Unknown "
                "MeshQualityType.");
            return;
        }
        _quality_tester->calculateQuality();
    }

    MeshLib::MeshQualityType const _type;
    MeshLib::Mesh const& _mesh;
    std::unique_ptr<MeshToolsLib::ElementQualityMetric> _quality_tester;
};

}  // namespace MeshToolsLib
