/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cassert>
#include <memory>

namespace BaseLib
{
class ConfigTree;
}

namespace MaterialLib
{
namespace Fracture
{
struct PermeabilityState
{
    virtual ~PermeabilityState() = default;
};

/**
 * Interface for fracture permeability models.
 */
class Permeability
{
public:
    virtual double permeability(PermeabilityState const* const state,
                                double const aperture0,
                                double const aperture_m) const = 0;

    virtual double dpermeability_daperture(PermeabilityState const* const state,
                                           double const aperture0,
                                           double const aperture_m) const = 0;

    virtual ~Permeability() = default;

    virtual std::unique_ptr<PermeabilityState> getNewState() const = 0;
};

/// A constant permeability model.
class ConstantPermeability final : public Permeability
{
public:
    ConstantPermeability(double const permeability)
        : _permeability(permeability)
    {
    }

private:
    double permeability(PermeabilityState const* const /*state*/,
                        double const /*aperture0*/,
                        double const /*aperture_m*/) const override
    {
        return _permeability;
    }

    double dpermeability_daperture(PermeabilityState const* const /*state*/,
                                   double const /*aperture0*/,
                                   double const /*aperture_m*/) const override
    {
        return 0;
    }

    std::unique_ptr<PermeabilityState> getNewState() const override
    {
        return nullptr;
    }

private:
    double const _permeability;
};

std::unique_ptr<ConstantPermeability> createConstantPermeability(
    BaseLib::ConfigTree const& config);

/// Hydraulic aperture equals the initial mechanical aperture s.t. permeability
/// is constant.
class ConstantHydraulicAperture final : public Permeability
{
    double permeability(PermeabilityState const* const /*state*/,
                        double const aperture0,
                        double const /*aperture_m*/) const override
    {
        return aperture0 * aperture0 / 12;
    }

    double dpermeability_daperture(PermeabilityState const* const /*state*/,
                                   double const /*aperture0*/,
                                   double const /*aperture_m*/) const override
    {
        return 0;
    }

    std::unique_ptr<PermeabilityState> getNewState() const override
    {
        return nullptr;
    }
};

std::unique_ptr<ConstantHydraulicAperture> createConstantHydraulicAperture(
    BaseLib::ConfigTree const& config);

/// Hydraulic aperture equals the mechanical aperture s.t. multiplication of the
/// permeability by the mechanical aperture yields the cubic law.
class CubicLaw final : public Permeability
{
public:
    explicit CubicLaw(double const permeability_threshold)
        : _permeability_threshold(permeability_threshold)
    {
    }

private:
    double permeability(PermeabilityState const* const /*state*/,
                        double const /*aperture0*/,
                        double const aperture_m) const override
    {
        return (aperture_m < _permeability_threshold)
                   ? 0
                   : aperture_m * aperture_m / 12;
    }

    double dpermeability_daperture(PermeabilityState const* const /*state*/,
                                   double const /*aperture0*/,
                                   double const aperture_m) const override
    {
        return (aperture_m < _permeability_threshold) ? 0 : aperture_m / 6;
    }

    std::unique_ptr<PermeabilityState> getNewState() const override
    {
        return nullptr;
    }

private:
    /// Zero permeability below the threshold value.
    double const _permeability_threshold;
};

std::unique_ptr<CubicLaw> createCubicLaw(BaseLib::ConfigTree const& config);

struct CubicLawAfterShearSlipState : public PermeabilityState
{
    explicit CubicLawAfterShearSlipState(bool const shear_slip_occured_)
        : shear_slip_occured(shear_slip_occured_)
    {
    }

    bool shear_slip_occured;
};

/// Hydraulic aperture ...
class CubicLawAfterShearSlip final : public Permeability
{
public:
    explicit CubicLawAfterShearSlip(double const initial_creation_aperture,
                                    double const minimum_permeability,
                                    double const permeability_threshold)
        : _initial_creation_aperture(initial_creation_aperture),
          _minimum_permeability(minimum_permeability),
          _permeability_threshold(permeability_threshold)
    {
    }

private:
    bool shearSlipOccured(PermeabilityState const* const state) const
    {
        assert(state != nullptr);
        return static_cast<CubicLawAfterShearSlipState const*>(state)
            ->shear_slip_occured;
    }

    double permeability(PermeabilityState const* const state,
                        double const /*aperture0*/,
                        double const aperture_m) const override
    {
        if (aperture_m < _permeability_threshold || !shearSlipOccured(state))
        {
            return _minimum_permeability;
        }

        double const aperture = aperture_m + _initial_creation_aperture;
        return aperture * aperture / 12;
    }

    double dpermeability_daperture(PermeabilityState const* const state,
                                   double const /*aperture0*/,
                                   double const aperture_m) const override
    {
        if (aperture_m < _permeability_threshold || !shearSlipOccured(state))
        {
            return 0;
        }

        double const aperture = aperture_m + _initial_creation_aperture;
        return aperture / 6;
    }

    std::unique_ptr<PermeabilityState> getNewState() const override
    {
        return std::make_unique<CubicLawAfterShearSlipState>(false);
    }

private:
    /// An induced initial creation aperture. As soon as shear- or tensile-
    /// failure occur, an induced creation aperture is added.
    double const _initial_creation_aperture;

    double const _minimum_permeability;

    /// Zero permeability below the threshold value.
    double const _permeability_threshold;
};

std::unique_ptr<CubicLawAfterShearSlip> createCubicLawAfterShearSlip(
    BaseLib::ConfigTree const& config);
}  // namespace Fracture
}  // namespace MaterialLib
