// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

/**
 * \file
 * \brief Interface for coupling OpenGeoSys with an external geochemical solver.
 *
 * Reactive transport in OGS is operator-split:
 *  - Transport: advection / diffusion of mobile (dissolved) components
 *    (PDE solve in OGS).
 *  - Chemistry: local reaction / speciation solved by an external
 *    chemistry backend (e.g. PHREEQC) for many independent local
 *    systems.
 *
 * Each reactive control volume is identified by \c chemical_system_id.
 * A \c chemical_system_id typically corresponds to one integration
 * point of an element (i.e. one local reaction volume). During the
 * chemistry step each \c chemical_system_id is treated as a well-mixed
 * batch reactor. There is no mass exchange between different
 * \c chemical_system_id in this step; spatial coupling happens only
 * in the transport solve.
 *
 * Per-timestep sequence (time increment dt):
 *
 *  1. OGS solves flow / transport and updates, for each control
 *     volume:
 *       - total dissolved component amounts \f$c_{T\alpha}\f$
 *         (the transported primary components),
 *       - pore pressure p,
 *       - temperature T (if a thermal process is active),
 *       - porosity / saturation.
 *
 *  2. setChemicalSystemConcrete(...):
 *     For each \c chemical_system_id, OGS provides the chemistry
 *     backend with the current local state at that point:
 *       - transported component totals \f$c_{T\alpha}\f$,
 *       - pore pressure,
 *       - temperature,
 *       - porosity and saturation,
 *       - reactive mineral amounts and reactive surface area,
 *       - current time t and time step size dt.
 *     This populates the backend's internal representation of that
 *     local chemical system (e.g. PHREEQC SOLUTION,
 *     EQUILIBRIUM_PHASES / KINETICS, SURFACE, EXCHANGE blocks).
 *
 *  3. executeSpeciationCalculation(dt):
 *     The backend advances geochemistry for all
 *     \c chemical_system_id over dt. Each local system is advanced as
 *     a closed batch reactor. The backend computes:
 *       - updated aqueous composition,
 *       - reaction/source terms \f$R_{\alpha}\f$ for each transported
 *         component,
 *       - updated mineral amounts.
 *
 *  4. updateVolumeFractionPostReaction(...),
 *     updatePorosityPostReaction(...):
 *     Reaction-induced precipitation / dissolution is pushed back to
 *     OGS:
 *       - solid volume fractions and mineral inventories are updated,
 *       - porosity is updated,
 *       - permeability / other medium properties can be updated
 *         based on the new porosity.
 *     These updated medium properties are then used in the next flow
 *     / transport step.
 *
 *  5. computeSecondaryVariable(...):
 *     Assemble derived output fields (e.g. pH, mineral volume
 *     fractions, surface loading) for writing to the mesh.
 *
 * Responsibilities of a concrete ChemicalSolverInterface
 * implementation (e.g. the PHREEQC-based implementation):
 *  - store and initialize per-\c chemical_system_id chemical state
 *    (component inventories, pH / pe, kinetic and equilibrium
 *    reactants, ion-exchange / surface sites, etc.);
 *  - accept per-\c chemical_system_id inputs from OGS before each
 *    reaction step (component totals, pressure, temperature,
 *    porosity, etc.);
 *  - run local reaction / speciation over dt;
 *  - provide reaction/source terms and updated mineral / porosity
 *    data back to OGS;
 *  - expose component totals and component names for the next
 *    transport solve.
 *
 * Thermodynamic state:
 *  Pressure and temperature are expected to be provided per
 *  \c chemical_system_id.
 */
#pragma once

#include <Eigen/Sparse>

#include "MaterialLib/MPL/VariableType.h"
#include "MathLib/LinAlg/GlobalLinearSolverType.h"
#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"
#include "MeshLib/Mesh.h"

namespace MaterialPropertyLib
{
class Medium;
}

namespace ParameterLib
{
class SpatialPosition;
}

namespace ChemistryLib
{
class ChemicalSolverInterface
{
public:
    ChemicalSolverInterface(MeshLib::Mesh const& mesh,
                            GlobalLinearSolver& linear_solver_)
        : _mesh(mesh), linear_solver(linear_solver_)
    {
        auto const* const bulk_element_ids = bulkElementIDs(_mesh);
        if (bulk_element_ids == nullptr)
        {
            OGS_FATAL(
                "The 'bulk_element_ids' property does not exist on the mesh "
                "{:s}.",
                _mesh.getName());
        }
        active_element_ids_.assign(bulk_element_ids->begin(),
                                   bulk_element_ids->end());
    }

    std::vector<std::size_t> const& activeElementIDs() const
    {
        return active_element_ids_;
    }

    virtual void initialize() {}

    virtual void initializeChemicalSystemConcrete(
        std::vector<double> const& /*concentrations*/,
        GlobalIndexType const& /*chemical_system_id*/,
        MaterialPropertyLib::Medium const& /*medium*/,
        ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/)
    {
    }

    virtual void setChemicalSystemConcrete(
        std::vector<double> const& /*concentrations*/,
        GlobalIndexType const& /*chemical_system_id*/,
        MaterialPropertyLib::Medium const* /*medium*/,
        MaterialPropertyLib::VariableArray const& /*vars*/,
        ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
        double const /*dt*/)
    {
    }

    virtual void setAqueousSolutionsPrevFromDumpFile() {}

    /**
     * Run PHREEQC for all local chemical systems for the current timestep.
     *
     * Uses the state previously provided by initializeChemicalSystemConcrete()
     * / setChemicalSystemConcrete() for each \c chemical_system_id. Each local
     * system is advanced as an isolated batch reactor over \f$\Delta t\f$:
     * no mass is exchanged between different systems inside PHREEQC.
     *
     * PHREEQC returns, for each system:
     *  - updated aqueous composition,
     *  - reaction / source terms \f$R_{\alpha}\f$ for each transported
     *    component,
     *  - updated mineral amounts and saturation state.
     *
     * After this call completes, these reacted values are available for
     * porosity update and for assembling the reaction/source term in the
     * transport equation.
     *
     * \param dt  Timestep size \f$\Delta t\f$ used for kinetic reactions.
     */
    virtual void executeSpeciationCalculation([[maybe_unused]] double const dt)
    {
    }

    virtual double getConcentration(
        int const /*component_id*/,
        GlobalIndexType const /*chemical_system_id*/) const
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    virtual std::vector<std::string> const getComponentList() const
    {
        return {};
    }

    virtual Eigen::SparseMatrix<double> const* getStoichiometricMatrix() const
    {
        return nullptr;
    }

    virtual double getKineticPrefactor(
        [[maybe_unused]] std::size_t reaction_id) const
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    /**
     * Apply mineral precipitation / dissolution to the solid fraction in OGS.
     *
     * Called after executeSpeciationCalculation(). For a given
     * \c chemical_system_id, this function updates the solid / mineral
     * inventory and related volume fractions based on the PHREEQC results
     * for the last timestep.
     *
     * Typical effect:
     *  - change in mineral volume fraction at this location,
     *  - updated solid composition available to mechanics / flow.
     */
    virtual void updateVolumeFractionPostReaction(
        GlobalIndexType const& /*chemical_system_id*/,
        MaterialPropertyLib::Medium const& /*medium*/,
        ParameterLib::SpatialPosition const& /*pos*/, double const /*porosity*/,
        double const /*t*/, double const /*dt*/)
    {
    }

    /**
     * Update porosity after chemical reactions.
     *
     * This applies porosity changes caused by precipitation / dissolution in
     * the local chemical system. The updated porosity can then be used to
     * update flow parameters (e.g. permeability) in the next timestep.
     */
    virtual void updatePorosityPostReaction(
        GlobalIndexType const& /*chemical_system_id*/,
        MaterialPropertyLib::Medium const& /*medium*/,
        double& /*porosity*/)
    {
    }

    virtual void computeSecondaryVariable(
        std::size_t const /*ele_id*/,
        std::vector<GlobalIndexType> const& /*chemical_system_indices*/)
    {
    }

    virtual ~ChemicalSolverInterface() = default;

public:
    std::vector<GlobalIndexType> chemical_system_index_map;
    MeshLib::Mesh const& _mesh;
    /// specify the linear solver used to solve the linearized reaction
    /// equation.
    GlobalLinearSolver& linear_solver;

private:
    std::vector<std::size_t> active_element_ids_;
};
}  // namespace ChemistryLib
