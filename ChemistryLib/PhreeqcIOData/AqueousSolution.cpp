/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "AqueousSolution.h"

#include <cmath>
#include <ostream>

#include "ChemistryLib/Common/ChargeBalance.h"
#include "MathLib/LinAlg/LinAlg.h"
#include "MeshLib/PropertyVector.h"

namespace ChemistryLib
{
namespace PhreeqcIOData
{
void AqueousSolution::print(std::ostream& os,
                            std::size_t const chemical_system_id) const
{
    os << "temp " << temperature << "\n";

    os << "pressure " << pressure << "\n";

    MathLib::LinAlg::setLocalAccessibleVector(*pH);
    switch (charge_balance)
    {
        case ChargeBalance::pH:
            os << "pH " << -std::log10((*pH)[chemical_system_id]) << " charge"
               << "\n";
            os << "pe " << (*pe)[chemical_system_id] << "\n";
            break;
        case ChargeBalance::pe:
            os << "pH " << -std::log10((*pH)[chemical_system_id]) << "\n";
            os << "pe " << (*pe)[chemical_system_id] << " charge"
               << "\n";
            break;
        case ChargeBalance::Unspecified:
            os << "pH " << -std::log10((*pH)[chemical_system_id]) << "\n";
            os << "pe " << (*pe)[chemical_system_id] << "\n";
            break;
    }

    os << "units mol/kgw\n";

    for (auto const& component : components)
    {
        MathLib::LinAlg::setLocalAccessibleVector(*component.amount);
        os << component.name << " " << (*component.amount)[chemical_system_id];
        component.chemical_formula.empty()
            ? os << "\n"
            : os << " as " << component.chemical_formula << "\n";
    }

    os << "\n\n";
}
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
