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

#include <limits>
#include <vector>

namespace ChemistryLib
{
namespace SelfContainedSolverData
{
struct ChemicalReaction
{
    explicit ChemicalReaction(std::vector<double> stoichiometric_vector_)
        : stoichiometric_vector(std::move(stoichiometric_vector_))
    {
    }

    virtual ~ChemicalReaction() = default;

    virtual double getKineticPrefactor() const
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    std::vector<double> stoichiometric_vector;
};

/**
 * \brief  First-order reaction
 * \details The rate law of the first-order reaction has the form of
 * \f[ Rate = k * c \f]
 * \f$ k \f$ is the first-order rate constant [1/s], and
 * \f$ c \f$ is the concentration [mol/m3]
 */
struct FirstOrderReaction : public ChemicalReaction
{
    FirstOrderReaction(std::vector<double> stoichiometric_vector_,
                       double first_order_rate_constant_)
        : ChemicalReaction(std::move(stoichiometric_vector_)),
          first_order_rate_constant(first_order_rate_constant_)
    {
    }

    double getKineticPrefactor() const override
    {
        return first_order_rate_constant;
    }

    /// the first order rate constant for first-order reaction.
private:
    double first_order_rate_constant;
};
}  // namespace SelfContainedSolverData
}  // namespace ChemistryLib
