// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

/**
 * \file
 * \brief Per-system aqueous state exchanged with PHREEQC.
 *
 * OpenGeoSys maps each reactive control volume to a local chemical system
 * identified externally by \c chemical_system_id. Each local chemical system
 * is treated as a closed, well-mixed batch reactor during the chemistry step:
 * no mass is exchanged with any other system while PHREEQC is running.
 *
 * This file defines two core data structures:
 *
 *  - Component:
 *    Represents one transported component (e.g. Na, Cl, Ca). For each
 *    \c chemical_system_id, \c Component::amount holds the total inventory
 *    of that component in that local chemical system. These totals
 *    \f$c_{T\alpha}\f$ are passed to PHREEQC as input and are updated after
 *    the chemistry step.
 *
 *  - AqueousSolution:
 *    Collects the input state needed to build the PHREEQC SOLUTION block
 *    for one local chemical system:
 *      - temperature \f$T\f$ [K],
 *      - pressure \f$p\f$ [Pa],
 *      - redox control (pe, \c fixing_pe, \c pe0),
 *      - charge-balance mode,
 *      - the per-component totals \f$c_{T\alpha}\f$.
 *
 *    After PHREEQC runs, AqueousSolution is also used to store the reacted
 *    output for that same \c chemical_system_id (updated pH, pe, etc.).
 *
 * Notes:
 *  - \c Component::amount, pH, and pe are stored per \c chemical_system_id.
 *    After each PHREEQC call, these per-system values are updated and then
 *    written back into OpenGeoSys for the next transport step.
 */
#pragma once

#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"
#include "Output.h"

namespace MeshLib
{
template <typename PROP_VAL_TYPE>
class PropertyVector;
}

namespace ChemistryLib
{
enum class ChargeBalance;

namespace PhreeqcIOData
{
struct Component
{
    explicit Component(std::string name_, std::string chemical_formula_)
        : name(std::move(name_)), chemical_formula(std::move(chemical_formula_))
    {
    }

    std::string const name;
    std::string const chemical_formula;
    std::unique_ptr<GlobalVector> amount;
    static const ItemType item_type = ItemType::Component;
};

struct AqueousSolution
{
    AqueousSolution(bool const fixing_pe_, double temperature_,
                    double pressure_, MeshLib::PropertyVector<double>* pe_,
                    double const pe0_, std::vector<Component>&& components_,
                    ChargeBalance charge_balance_)
        : fixing_pe(fixing_pe_),
          temperature(temperature_),
          pressure(pressure_),
          pe(pe_),
          pe0(pe0_),
          components(std::move(components_)),
          charge_balance(charge_balance_)
    {
    }

    void print(std::ostream& os, std::size_t const chemical_system_id) const;

    /// When this option is enabled, the pe value will be fixed over time by
    /// adding or removing atmospheric oxygen.
    bool const fixing_pe;
    double const temperature;
    double const pressure;
    std::unique_ptr<GlobalVector> pH;
    MeshLib::PropertyVector<double>* pe;
    double const pe0;
    std::vector<Component> components;
    ChargeBalance const charge_balance;
};
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
