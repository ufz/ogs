/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file:   ComponentProperties.h
 */

#pragma once

#include <memory>

#include <ProcessLib/Parameter/Parameter.h>

namespace MaterialLib
{
namespace Component
{
class ComponentProperties final
{
public:
    ComponentProperties(
        std::string const name_,
        ProcessLib::Parameter<double> const& molecular_diffusion_,
        ProcessLib::Parameter<double> const& solute_dispersivity_longitudinal_,
        ProcessLib::Parameter<double> const& solute_dispersivity_transverse_)
        : _name(name_),
          _molecular_diffusion(molecular_diffusion_),
          _solute_dispersivity_longitudinal(solute_dispersivity_longitudinal_),
          _solute_dispersivity_transverse(solute_dispersivity_transverse_)
    {
    }

    double getMolecularDiffusion(double const t,
                                 ProcessLib::SpatialPosition const& x) const;

    double getLongitudinalDispersivity(
        double const t, ProcessLib::SpatialPosition const& x) const;

    double getTransversalDispersivity(
        double const t, ProcessLib::SpatialPosition const& x) const;

private:
    std::string const _name;
    ProcessLib::Parameter<double> const& _molecular_diffusion;
    ProcessLib::Parameter<double> const& _solute_dispersivity_longitudinal;
    ProcessLib::Parameter<double> const& _solute_dispersivity_transverse;
};
}  // namespace Component
}  // namespace MaterialLib
