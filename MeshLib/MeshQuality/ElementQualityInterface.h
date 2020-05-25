/**
 * \file
 * \author Karsten Rink
 * \date   2015-03-24
 * \brief  Definition of the ElementQualityInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "BaseLib/Histogram.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/MeshQuality/ElementQualityMetric.h"
#include "MeshLib/MeshQuality/EdgeRatioMetric.h"
#include "MeshLib/MeshQuality/ElementSizeMetric.h"
#include "MeshLib/MeshQuality/SizeDifferenceMetric.h"
#include "MeshLib/MeshQuality/AngleSkewMetric.h"
#include "MeshLib/MeshQuality/RadiusEdgeRatioMetric.h"

namespace MeshLib
{

/**
 * Interface class for handling mesh element quality metrics
 */
class ElementQualityInterface
{
public:
    /// Constructor
    ElementQualityInterface(MeshLib::Mesh const& mesh, MeshQualityType t)
    : type_(t), mesh_(mesh), quality_tester_(nullptr)
    {
        calculateElementQuality(mesh_, type_);
    }

    /// Destructor
    ~ElementQualityInterface()
    {
        delete quality_tester_;
    }

    /// Returns the vector containing a quality measure for each element.
    std::vector<double> const getQualityVector() const
    {
        if (quality_tester_)
            return quality_tester_->getElementQuality();

        std::vector<double> empty_quality_vec(0);
        return empty_quality_vec;
    }

    /// Returns a histogram of the quality vector separated into the given number of bins.
    /// If no number of bins is specified, one will be calculated based on the Sturges criterium.
    BaseLib::Histogram<double> getHistogram(std::size_t n_bins = 0) const
    {
        if (quality_tester_)
            return quality_tester_->getHistogram(static_cast<std::size_t>(n_bins));

        return BaseLib::Histogram<double>{{}};
    }

    /// Writes a histogram of the quality vector to a specified file.
    int writeHistogram(std::string const& file_name, std::size_t n_bins = 0) const
    {
        if (quality_tester_ == nullptr)
            return 1;

        BaseLib::Histogram<double> const histogram (quality_tester_->getHistogram(n_bins));
        histogram.write(file_name, mesh_.getName(), MeshQualityType2String(type_));
        return 0;
    }

private:
    /// Calculates the quality of each mesh element based on the specified metric
    void calculateElementQuality(MeshLib::Mesh const& mesh, MeshQualityType t)
    {
        if (t == MeshQualityType::EDGERATIO)
            quality_tester_ = new MeshLib::EdgeRatioMetric(mesh);
        else if (t == MeshQualityType::ELEMENTSIZE)
            quality_tester_ = new MeshLib::ElementSizeMetric(mesh);
        else if (t == MeshQualityType::SIZEDIFFERENCE)
            quality_tester_ = new MeshLib::SizeDifferenceMetric(mesh);
        else if (t == MeshQualityType::EQUIANGLESKEW)
            quality_tester_ = new MeshLib::AngleSkewMetric(mesh);
        else if (t == MeshQualityType::RADIUSEDGERATIO)
            quality_tester_ = new MeshLib::RadiusEdgeRatioMetric(mesh);
        else
        {
            ERR("ElementQualityInterface::calculateElementQuality(): Unknown MeshQualityType.");
            return;
        }
        quality_tester_->calculateQuality();
    }

    MeshQualityType const type_;
    MeshLib::Mesh const& mesh_;
    MeshLib::ElementQualityMetric* quality_tester_;
};

}  // namespace MeshLib
