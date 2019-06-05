/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <fstream>

#include "AqueousSolution.h"

namespace ChemistryLib
{
std::ofstream& operator<<(std::ofstream& out,
                          AqueousSolution const& aqueous_solution)
{
    out << "temp " << aqueous_solution.temperature << "\n";

    out << "pressure " << aqueous_solution.pressure << "\n";

    switch (aqueous_solution.means_of_adjusting_charge)
    {
        case MeansOfAdjustingCharge::pH:
            out << "pH " << aqueous_solution.pH << " charge"
                << "\n";
            out << "pe " << aqueous_solution.pe << "\n";
            break;
        case MeansOfAdjustingCharge::pe:
            out << "pH " << aqueous_solution.pH << "\n";
            out << "pe " << aqueous_solution.pe << " charge"
                << "\n";
            break;
        case MeansOfAdjustingCharge::Unspecified:
            out << "pH " << aqueous_solution.pH << "\n";
            out << "pe " << aqueous_solution.pe << "\n";
            break;
    }

    out << "units mol/kgw\n";

    for (auto const& component : aqueous_solution.components)
    {
        out << component.name << " " << component.amount << "\n";
    }

    return out;
}
}  // namespace ChemistryLib
