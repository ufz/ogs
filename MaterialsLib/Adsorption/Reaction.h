/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include<memory>

#include "BaseLib/ConfigTreeNew.h"

namespace Adsorption
{

class Reaction
{
public:
	static std::unique_ptr<Reaction> newInstance(BaseLib::ConfigTreeNew const& rsys);

	virtual double get_enthalpy(const double p_Ads, const double T_Ads, const double M_Ads) const = 0;
	virtual double get_reaction_rate(const double p_Ads, const double T_Ads,
									 const double M_Ads, const double loading) const = 0;

	// TODO: get rid of
	virtual double get_equilibrium_loading(const double, const double, const double) const {
		return -1.0;
	}

    virtual ~Reaction() = default;
};

}
