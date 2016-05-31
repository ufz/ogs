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
#include "ProcessLib/Process.h" // TODO move findParam somewhere else

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
    IdealGasBinaryMixtureDensity(
        BaseLib::ConfigTree const& config,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters)
        : _M0(ProcessLib::findParameter<double>(config, "molar_mass_0", parameters)),
          _M1(ProcessLib::findParameter<double>(config, "molar_mass_1", parameters))
    {
        config.checkConfParam("type", "IdealGasBinaryMixtureDensity");
    }

    double getValue(const double t, const double* const x,
                    const GlobalIndexType node,
                    const MeshLib::Element& element,
                    const std::size_t integration_point, const double& p,
                    const double& T, const double& x_m) const
    {
        const double M0 = _M0(t, x, node, element, integration_point);
        const double M1 = _M1(t, x, node, element, integration_point);
        const double xn = M0 * x_m / (M0 * x_m + M1 * (1.0 - x_m));

        return p / (GAS_CONST * T) * (M1 * xn + M0 * (1.0 - xn));
    }

    double getDerivative(const unsigned derivative_in_direction_of_argument,
                         const double /*t*/, const double* const /*x*/,
                         const GlobalIndexType /*node*/,
                         const MeshLib::Element& /*element*/,
                         const std::size_t /*integration_point*/,
                         const double& p, const double& T,
                         const double& x) const
    {
        (void)p;
        (void)T;
        (void)x;
        (void)derivative_in_direction_of_argument;
        ERR("TODO implement.");
        std::abort();
        return 0.0;
    }

private:
    ProcessLib::Parameter<double> const& _M0;
    ProcessLib::Parameter<double> const& _M1;
};

OGS_DEFINE_CONSTITUTIVE_RELATION_BUILDER(IdealGasBinaryMixtureDensity, double,
                                         double, double, double)

}  // namespace ConstitutiveRelation
}  // namespace ProcessLib
