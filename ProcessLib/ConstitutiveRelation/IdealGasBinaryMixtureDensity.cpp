/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "IdealGasBinaryMixtureDensity.h"

#include "BaseLib/ConfigTree.h"

namespace
{
const double GAS_CONST = 8.3144621;
}

namespace ProcessLib
{
namespace ConstitutiveRelation
{
class IdealGasBinaryMixtureDensity
    : public ConstitutiveRelation<double, double, double, double>
{
public:
    IdealGasBinaryMixtureDensity(BaseLib::ConfigTree const& config)
        : _M0(config.getConfParam<double>("M0")),
          _M1(config.getConfParam<double>("M1"))
    {
    }

    double getValue(const double /*t*/, const double* const /*x*/,
                    const GlobalIndexType /*node*/,
                    const MeshLib::Element& /*element*/,
                    const std::size_t /*integration_point*/, const double& p,
                    const double& T, const double& x) const
    {
        const double xn = _M0 * x / (_M0 * x + _M1 * (1.0 - x));

        return p / (GAS_CONST * T) * (_M1 * xn + _M0 * (1.0 - xn));
    }

    double getDerivative(const unsigned derivative_in_direction_of_argument,
                         const double /*t*/, const double* const /*x*/,
                         const GlobalIndexType /*node*/,
                         const MeshLib::Element& /*element*/,
                         const std::size_t /*integration_point*/,
                         const double& p, const double& T,
                         const double& x) const
    {
        (void) p; (void) T; (void) x; (void) derivative_in_direction_of_argument;
        ERR("TODO implement.");
        std::abort();
        return 0.0;
    }

private:
    const double _M0;
    const double _M1;
};

std::unique_ptr<ConstitutiveRelation<double, double, double, double>>
IdealGasBinaryMixtureDensityBuilder::createConstitutiveRelation(
    BaseLib::ConfigTree const& config)
{
    return std::unique_ptr<
        ConstitutiveRelation<double, double, double, double>>(
            new IdealGasBinaryMixtureDensity(config));
}

}  // namespace ConstitutiveRelation
}  // namespace ProcessLib
