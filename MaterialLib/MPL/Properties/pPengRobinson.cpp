/**
 * \author Norbert Grunwald
 * \date   18.09.2017
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "pPengRobinson.h"
#include <algorithm>
#include <cmath>
#include "../../MathLib/poly34.h"
#include "../mpComponent.h"
#include "../mpPhase.h"
#include "pUniversalConstants.h"

namespace MaterialPropertyLib
{
/// This constructor throws an error, since the property is not
/// implemented on the medium scale.
PengRobinson::PengRobinson(Medium* m) : _phase(0), _component(0)
{
    notImplemented("PengRobinson", "Medium");
}
/// Constructor for the phase version (binary mixture) of
/// the Peng-Robinson EOS
PengRobinson::PengRobinson(Phase* p) : _phase(p), _component(0){};

/// Constructor for the component version (pure substances) of
/// the Peng-Robinson EOS
PengRobinson::PengRobinson(Component* c) : _phase(0), _component(c){};

/**
 */
PropertyDataType PengRobinson::value(VariableArray const& v)
{
    if (isUpdated())
        return _value;

    double M(0);
    double k_ij(0);
    double am(0), bm(0);

    const double temperature = getScalar(v[T]);

    // A temporary component vector holds all the properties from
    // either the (up to two) phase components (if used as phase
    // property), or of the component (if used as component
    // property).
    std::vector<Component*> temp_components;

    if (_phase)
    {
        // Implementation as phase property: EOS for binary mixtures or
        // for single fluids. If more than two components are defined, only
        // the first two components act as primary components, all others
        // will be neglected.
        const auto nComponents = std::min(
            std::max(1, static_cast<int>(_phase->numberOfComponents())), 2);
        // The number of relevant components is now narrowed down to
        // one or two components.
        for (int c = 0; c < nComponents; ++c)
            temp_components.push_back(&_phase->component(c));

        M = getScalar(_phase->property(molar_mass));
        k_ij = getScalar(_phase->property(binary_interaction_coefficient));
    }
    else
    {
        // Implementation as component property. The component properties
        // are simply copied to the temporary component vector.
        temp_components.push_back(_component);
        M = getScalar(_component->property(molar_mass));
        k_ij = 0;
    }

    // from now on, the implementation is identical for both phase and
    // component property...

    const std::size_t nComponents = temp_components.size();

    assert(nComponents == 1 || nComponents == 2);

    std::array<double, 2> sqrt_a, b, x_n;

    // In order to solve the EOS, we need to obtain two characteristic
    // parameters a and b. In case of a binary mixture, these parameters
    // are computed for each compound separately and later combined by
    // some mixing rule.

    for (std::size_t c = 0; c < nComponents; ++c)
    {
        const double T_crit =
            getScalar(temp_components[c]->property(critical_temperature));
        const double p_crit =
            getScalar(temp_components[c]->property(critical_pressure));
        const double omega =
            getScalar(temp_components[c]->property(acentric_factor));

        const double alpha_Tr = alpha(temperature, T_crit, omega);
        const double a_Tc = cohesionPressure(T_crit, p_crit);

        sqrt_a[c] = std::sqrt(a_Tc * alpha_Tr);
        b[c] = coVolume(T_crit, p_crit);
        x_n[c] = getScalar(temp_components[c]->property(mole_fraction));
    }

    for (std::size_t i = 0; i < nComponents; ++i)
        for (std::size_t j = 0; j < nComponents; ++j)
        {
            am += x_n[i] * x_n[j] * (1 - k_ij) * sqrt_a[i] * sqrt_a[j];
            bm += x_n[i] * b[i];
        }

    const double pressure = getScalar(v[p_GR]);
    const double RT = gasConstant * temperature;

    const double A = am * pressure / RT / RT;
    const double B = bm * pressure / RT;

    // Peng-Robinson EOS in terms of compressibility factor Z:
    // Z^3 + p*Z^2 + q*Z + r = 0
    // where p, q, r are defined by
    const double p = B - 1.;
    const double q = A - 3. * B * B - 2. * B;
    const double r = B * B * B + B * B - A * B;

    // The PR-EOS has at most three real roots.
    std::array<double, 3> roots;
    const int numberRoots = SolveP3(roots.data(), p, q, r);

    // This implementation only works for the single phase region (i.e. the
    // free spaces of the phase diagram, not to be confused with single- or
    // multiphase flow).

    const double Z =
        (numberRoots > 1 ? *std::max_element(std::begin(roots), std::end(roots))
                         : roots[0]);

    _value = M * pressure / Z / gasConstant / temperature;

    isUpdated(true);

    return _value;
}
double kappa(const double omega)
{
    return 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
}

double alpha(const double T, const double T_c, const double omega)
{
    const double sqrt_alpha = 1 + kappa(omega) * (1 - std::sqrt(T / T_c));
    return sqrt_alpha * sqrt_alpha;
}

double cohesionPressure(const double Tc, const double pc)
{
    const double R = gasConstant;
    return 0.45723553 * R * R * Tc * Tc / pc;
}

double coVolume(const double Tc, const double pc)
{
    const double R = gasConstant;
    return 0.077796074 * R * Tc / pc;
}

}  // MaterialPropertyLib
