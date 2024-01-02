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

#include "GeoLib/Raster.h"
#include "Parameter.h"

namespace ParameterLib
{
struct RasterParameter final : public Parameter<double>
{
    explicit RasterParameter(std::string const& name_,
                             GeoLib::NamedRaster const& named_raster)
        : Parameter<double>(name_), _named_raster(named_raster)
    {
    }

    bool isTimeDependent() const override { return false; }

    int getNumberOfGlobalComponents() const override { return 1; }

    std::vector<double> operator()(double const /*t*/,
                                   SpatialPosition const& pos) const override
    {
        auto const& coordinates = pos.getCoordinates();
        if (!coordinates)
        {
            OGS_FATAL("RasterParameter::operator(): couldn't get coordinates.");
        }
        auto const value = _named_raster.raster->getValueAtPoint(*coordinates);
        return {value};
    }

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
    getNodalValuesOnElement(MeshLib::Element const& element,
                            double const t) const override
    {
        auto const n_nodes = element.getNumberOfNodes();
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> result(
            n_nodes, getNumberOfGlobalComponents());

        SpatialPosition position;
        auto const nodes = element.getNodes();
        for (unsigned i = 0; i < n_nodes; ++i)
        {
            position.setCoordinates(*(nodes[i]));
            auto const& values = this->operator()(t, position);
            result.row(i) =
                Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1> const>(
                    values.data(), values.size());
        }

        return result;
    }

private:
    GeoLib::NamedRaster const& _named_raster;
};

std::unique_ptr<ParameterBase> createRasterParameter(
    std::string const& name, BaseLib::ConfigTree const& config,
    std::vector<GeoLib::NamedRaster> const& named_rasters);

}  // namespace ParameterLib
