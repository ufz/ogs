// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <IPhreeqc.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef USE_PETSC
#include "BaseLib/MPI.h"
#endif

#include "ChemistryLib/Common/ChargeBalance.h"
#include "ChemistryLib/PhreeqcIO.h"
#include "ChemistryLib/PhreeqcIOData/AqueousSolution.h"
#include "ChemistryLib/PhreeqcIOData/ChemicalSystem.h"
#include "ChemistryLib/PhreeqcIOData/CreateOutput.h"
#include "ChemistryLib/PhreeqcIOData/Dump.h"
#include "ChemistryLib/PhreeqcIOData/EquilibriumReactant.h"
#include "ChemistryLib/PhreeqcIOData/Exchange.h"
#include "ChemistryLib/PhreeqcIOData/KineticReactant.h"
#include "ChemistryLib/PhreeqcIOData/Knobs.h"
#include "ChemistryLib/PhreeqcIOData/Output.h"
#include "ChemistryLib/PhreeqcIOData/PhreeqcInstancePool.h"
#include "ChemistryLib/PhreeqcIOData/ReactionRate.h"
#include "ChemistryLib/PhreeqcIOData/UserPunch.h"
#include "InfoLib/TestInfo.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/VariableType.h"
#include "MathLib/LinAlg/Eigen/EigenOption.h"
#include "MathLib/LinAlg/GlobalLinearSolverType.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEnums.h"
#include "MeshLib/PropertyVector.h"
#include "MeshToolsLib/MeshGenerators/MeshGenerator.h"
#include "ParameterLib/SpatialPosition.h"
#include "Tests/MaterialLib/TestMPL.h"

using namespace ChemistryLib;
using namespace ChemistryLib::PhreeqcIOData;

namespace
{
std::string const phreeqc_database =
    TestInfoLib::TestInfo::data_path +
    "/Parabolic/ComponentTransport/ReactiveTransport/CationExchange/"
    "phreeqc.dat";

/// Return a default-constructed linear solver of the build's GlobalLinearSolver
/// type.  PhreeqcIO stores the reference but does not use it for the speciation
/// path exercised by these tests.
std::unique_ptr<GlobalLinearSolver> makeDummyLinearSolver()
{
#if defined(USE_LIS) || defined(USE_PETSC)
    return std::make_unique<GlobalLinearSolver>("test", "");
#else
    return std::make_unique<GlobalLinearSolver>("test", MathLib::EigenOption{});
#endif
}

/// Test driver that builds a PhreeqcIO over a tiny in-memory line mesh and
/// exposes a small API for setting per-cell inputs and running the production
/// PhreeqcIO::executeSpeciationCalculation path (parallel pool path is taken
/// because stream mode is enabled).
struct PhreeqcIOTestDriver
{
    struct ComponentSpec
    {
        std::string name;
        std::string chemical_formula;
    };
    struct EquilibriumSpec
    {
        std::string name;
        double saturation_index;
        double initial_molality;
    };
    struct ExchangerSpec
    {
        std::string name;
        double initial_molality;
    };
    struct KineticSpec
    {
        std::string name;
        std::string chemical_formula;
        double initial_molality;
        std::vector<std::string> rate_script;
    };

    explicit PhreeqcIOTestDriver(int num_systems_, int num_threads_ = 1)
        : num_systems(num_systems_), num_threads(num_threads_)
    {
    }

    PhreeqcIOTestDriver& addComponent(std::string name,
                                      std::string chemical_formula = "")
    {
        components_spec.push_back(
            {std::move(name), std::move(chemical_formula)});
        return *this;
    }
    PhreeqcIOTestDriver& setChargeBalance(ChargeBalance cb)
    {
        charge_balance = cb;
        return *this;
    }
    PhreeqcIOTestDriver& addEquilibriumPhase(std::string name,
                                             double saturation_index,
                                             double initial_molality)
    {
        equilibrium_spec.push_back(
            {std::move(name), saturation_index, initial_molality});
        return *this;
    }
    PhreeqcIOTestDriver& addExchanger(std::string name, double initial_molality)
    {
        exchanger_spec.push_back({std::move(name), initial_molality});
        return *this;
    }
    PhreeqcIOTestDriver& addKineticPhase(std::string name,
                                         std::string chemical_formula,
                                         double initial_molality,
                                         std::vector<std::string>
                                             rate_script)
    {
        kinetic_spec.push_back({std::move(name), std::move(chemical_formula),
                                initial_molality, std::move(rate_script)});
        return *this;
    }

    /// Build the mesh, ChemicalSystem, Output, Dump, Knobs, and PhreeqcIO.
    void build()
    {
        mesh.reset(MeshToolsLib::MeshGenerator::generateLineMesh(
            static_cast<unsigned>(num_systems), 1.0));

        // ChemicalSolverInterface requires a bulk_element_ids cell property.
        auto* bulk_ids =
            mesh->getProperties().createNewPropertyVector<std::size_t>(
                "bulk_element_ids", MeshLib::MeshItemType::Cell, 1);
        bulk_ids->resize(num_systems);
        std::iota(bulk_ids->begin(), bulk_ids->end(), 0);

        auto* pe = mesh->getProperties().createNewPropertyVector<double>(
            "pe", MeshLib::MeshItemType::Cell, 1);
        pe->resize(num_systems, 4.0);

        solver = makeDummyLinearSolver();

        std::vector<Component> aq_components;
        aq_components.reserve(components_spec.size());
        for (auto const& c : components_spec)
        {
            aq_components.emplace_back(c.name, c.chemical_formula);
        }

        auto aqueous = std::make_unique<AqueousSolution>(
            /*fixing_pe=*/false,
            // PHREEQC expects temperature in degrees Celsius and pressure in
            // atm (see AqueousSolution::print).  Matches the production .prj.
            /*temperature=*/25.0,
            /*pressure=*/1.0, pe,
            /*pe0=*/4.0, std::move(aq_components), charge_balance);

        std::vector<EquilibriumReactant> equilibrium_reactants;
        equilibrium_reactants.reserve(equilibrium_spec.size());
        for (auto const& e : equilibrium_spec)
        {
            auto* mol = createCellPropertyVector("eq_" + e.name + "_molality",
                                                 e.initial_molality);
            auto* mol_prev = createCellPropertyVector(
                "eq_" + e.name + "_molality_prev", e.initial_molality);
            auto* vf =
                createCellPropertyVector("eq_" + e.name + "_volume_frac", 0.0);
            auto* vf_prev = createCellPropertyVector(
                "eq_" + e.name + "_volume_frac_prev", 0.0);
            auto* mp_mol = createCellPropertyVector(
                "eq_" + e.name + "_mesh_prop_molality", e.initial_molality);
            equilibrium_reactants.emplace_back(
                e.name, mol, mol_prev, vf, vf_prev, mp_mol, e.saturation_index,
                /*reaction_irreversibility=*/"");
        }

        std::vector<ExchangeSite> exchangers;
        exchangers.reserve(exchanger_spec.size());
        for (auto const& x : exchanger_spec)
        {
            auto* mol = createCellPropertyVector("ex_" + x.name + "_molality",
                                                 x.initial_molality);
            exchangers.emplace_back(x.name, mol);
        }

        std::vector<KineticReactant> kinetic_reactants;
        std::vector<ReactionRate> reaction_rates;
        kinetic_reactants.reserve(kinetic_spec.size());
        reaction_rates.reserve(kinetic_spec.size());
        for (auto const& k : kinetic_spec)
        {
            auto* mol = createCellPropertyVector("kin_" + k.name + "_molality",
                                                 k.initial_molality);
            auto* mol_prev = createCellPropertyVector(
                "kin_" + k.name + "_molality_prev", k.initial_molality);
            auto* vf =
                createCellPropertyVector("kin_" + k.name + "_volume_frac", 0.0);
            auto* vf_prev = createCellPropertyVector(
                "kin_" + k.name + "_volume_frac_prev", 0.0);
            auto* mp_mol = createCellPropertyVector(
                "kin_" + k.name + "_mesh_prop_molality", k.initial_molality);
            kinetic_reactants.emplace_back(k.name, k.chemical_formula, mol,
                                           mol_prev, vf, vf_prev, mp_mol,
                                           /*parameters=*/std::vector<double>{},
                                           /*fix_amount=*/false);
            reaction_rates.emplace_back(
                k.name, std::vector<std::string>(k.rate_script));
        }

        auto chemical_system = std::make_unique<ChemicalSystem>(
            std::move(aqueous), std::move(kinetic_reactants),
            std::move(equilibrium_reactants), std::move(exchangers),
            std::vector<
                std::variant<DensityBasedSurfaceSite, MoleBasedSurfaceSite>>{});

        std::string const project_file_name = "TestPhreeqcParallel";
        auto output = createOutput(*chemical_system,
                                   /*user_punch=*/nullptr,
                                   /*use_high_precision=*/true,
                                   project_file_name);

        // Match production: Dump is wired only when exchangers or surface
        // sites exist (see CreateChemicalSolverInterface.cpp). This driver
        // does not model surface sites, so the exchanger flag alone gates it.
        auto dump = exchanger_spec.empty()
                        ? nullptr
                        : std::make_unique<Dump>(project_file_name);

        Knobs knobs{/*max_iterations=*/500,
                    /*relative_convergence_tolerance=*/1e-12,
                    /*tolerance=*/1e-15,
                    /*step_size=*/5,
                    /*scaling=*/false};

        phreeqc_io = std::make_unique<PhreeqcIO>(
            *mesh, *solver, project_file_name, std::string(phreeqc_database),
            std::move(chemical_system), std::move(reaction_rates),
            /*user_punch=*/nullptr, std::move(output), std::move(dump),
            std::move(knobs),
            /*use_stream_mode=*/true, num_threads,
            /*concentration_warning_threshold=*/-1e-12);

        // ChemicalSolverInterface::chemical_system_index_map is normally
        // populated by the transport process before initialize() runs; supply
        // a 1:1 mapping (cell i -> chemical_system_id i) so that
        // PhreeqcIO::initialize sees the right system count.
        phreeqc_io->chemical_system_index_map.resize(num_systems);
        std::iota(phreeqc_io->chemical_system_index_map.begin(),
                  phreeqc_io->chemical_system_index_map.end(), 0);

        phreeqc_io->initialize();

        // Build a MaterialPropertyLib::Medium that mirrors the chemical
        // system's reactant / exchanger inventory.  The constants
        // density = porosity = molar_volume = 1 collapse the production
        // molality = volume_fraction / density / porosity / molar_volume
        // back to molality = volume_fraction = initial_molality, so the
        // numerical outcome matches what the driver previously set by hand.
        medium = buildTestMedium();

        // Drive the production initialization path: populates
        // (*reactant.molality)[id] and (*reactant.volume_fraction)[id] from
        // the Medium, and (*exchanger.molality)[id] from the Medium's
        // <molality> property. Component amounts are placeholders here (the
        // real state is set later via setChemicalSystemConcrete), but the last
        // entry is the H+ activity 10^-pH: it must be a physically valid
        // positive value (pH 7 -> 1e-7), matching the production init path,
        // which reads it from the pH initial condition. A zero would mean
        // pH = +inf and trips the H+ guard in setAqueousSolution.
        std::vector<double> initial_concentrations(components_spec.size() + 1,
                                                   0.0);
        initial_concentrations.back() = std::pow(10.0, -7.0);
        for (int id = 0; id < num_systems; ++id)
        {
            ParameterLib::SpatialPosition pos;
            pos.setElementID(id);
            phreeqc_io->initializeChemicalSystemConcrete(
                initial_concentrations, id, *medium, pos, /*t=*/0.0);
        }
    }

    /// Drive the production entry point: pack `component_amounts` (zero for
    /// any missing transported component) plus the H+ activity 10^-pH into
    /// the `concentrations` vector and hand it to
    /// PhreeqcIO::setChemicalSystemConcrete, which calls setAqueousSolution
    /// and setReactantMolality for us.
    void setInitialState(int chemical_system_id,
                         std::map<std::string, double> const& component_amounts,
                         double pH, double dt = 1.0)
    {
        std::vector<double> concentrations(components_spec.size() + 1, 0.0);
        for (std::size_t i = 0; i < components_spec.size(); ++i)
        {
            auto const it = component_amounts.find(components_spec[i].name);
            if (it != component_amounts.end())
            {
                concentrations[i] = it->second;
            }
        }
        concentrations.back() = std::pow(10.0, -pH);

        MaterialPropertyLib::VariableArray vars;
        vars.porosity = 1.0;  // matches the Medium-level <porosity> below.
        ParameterLib::SpatialPosition pos;
        pos.setElementID(chemical_system_id);

        phreeqc_io->setChemicalSystemConcrete(concentrations,
                                              chemical_system_id, medium.get(),
                                              vars, pos, /*t=*/0.0, dt);
    }

    void run(double dt) { phreeqc_io->executeSpeciationCalculation(dt); }
    /// Read H+ activity through the production getter and convert it to pH.
    double pH(int chemical_system_id) const
    {
        return -std::log10(phreeqc_io->getConcentration(
            // H+ concentration comes after all chemical components.
            static_cast<int>(components_spec.size()), chemical_system_id));
    }
    double getConcentration(int component_id, int chemical_system_id) const
    {
        return phreeqc_io->getConcentration(component_id, chemical_system_id);
    }

    int const num_systems;
    int const num_threads;
    std::vector<ComponentSpec> components_spec;
    ChargeBalance charge_balance = ChargeBalance::Unspecified;
    std::vector<EquilibriumSpec> equilibrium_spec;
    std::vector<ExchangerSpec> exchanger_spec;
    std::vector<KineticSpec> kinetic_spec;

    std::unique_ptr<MeshLib::Mesh> mesh;
    std::unique_ptr<GlobalLinearSolver> solver;
    std::unique_ptr<PhreeqcIO> phreeqc_io;
    std::unique_ptr<MaterialPropertyLib::Medium> medium;

private:
    MeshLib::PropertyVector<double>* createCellPropertyVector(
        std::string const& name, double const initial_value)
    {
        auto* pv = mesh->getProperties().createNewPropertyVector<double>(
            name, MeshLib::MeshItemType::Cell, 1);
        pv->resize(num_systems, initial_value);
        return pv;
    }

    /// Emit an XML Medium with the constituents PhreeqcIO will look up.
    /// Each equilibrium / kinetic reactant becomes a Solid constituent
    /// carrying volume_fraction (= initial_molality) and molar_volume (= 1);
    /// each exchanger becomes a Solid constituent carrying a Constant
    /// molality (= initial_molality) consumed by initializeSiteMolality.
    std::unique_ptr<MaterialPropertyLib::Medium> buildTestMedium() const
    {
        std::stringstream xml;
        xml << "<medium>\n"
               "  <phases>\n"
               "    <phase>\n"
               "      <type>AqueousLiquid</type>\n"
               "      <properties>\n"
               "        <property>\n"
               "          <name>density</name>\n"
               "          <type>Constant</type>\n"
               "          <value>1</value>\n"
               "        </property>\n"
               "      </properties>\n"
               "    </phase>\n"
               "    <phase>\n"
               "      <type>Solid</type>\n"
               "      <components>\n";
        for (auto const& e : equilibrium_spec)
        {
            xml << "        <component>\n"
                << "          <name>" << e.name << "</name>\n"
                << "          <properties>\n"
                << "            <property><name>volume_fraction</name>"
                << "<type>Constant</type><value>" << e.initial_molality
                << "</value></property>\n"
                << "            <property><name>molar_volume</name>"
                << "<type>Constant</type><value>1</value></property>\n"
                << "          </properties>\n"
                << "        </component>\n";
        }
        for (auto const& k : kinetic_spec)
        {
            xml << "        <component>\n"
                << "          <name>" << k.name << "</name>\n"
                << "          <properties>\n"
                << "            <property><name>volume_fraction</name>"
                << "<type>Constant</type><value>" << k.initial_molality
                << "</value></property>\n"
                << "            <property><name>molar_volume</name>"
                << "<type>Constant</type><value>1</value></property>\n"
                << "          </properties>\n"
                << "        </component>\n";
        }
        for (auto const& x : exchanger_spec)
        {
            xml << "        <component>\n"
                << "          <name>" << x.name << "</name>\n"
                << "          <properties>\n"
                << "            <property><name>molality</name>"
                << "<type>Constant</type><value>" << x.initial_molality
                << "</value></property>\n"
                << "          </properties>\n"
                << "        </component>\n";
        }
        xml << "      </components>\n"
               "      <properties/>\n"
               "    </phase>\n"
               "  </phases>\n"
               "  <properties>\n"
               "    <property>\n"
               "      <name>porosity</name>\n"
               "      <type>Constant</type>\n"
               "      <value>1</value>\n"
               "    </property>\n"
               "  </properties>\n"
               "</medium>\n";
        return Tests::createTestMaterial(xml.str(), /*geometry_dimension=*/1);
    }
};

/// Build a driver with a dilute NaCl solution and a no-op Halite equilibrium
/// phase.  Halite (NaCl) at SI=-50 with mol=0 never precipitates, so it leaves
/// the chemistry untouched while still emitting an EQUILIBRIUM_PHASES block and
/// the "reacted" SELECTED_OUTPUT row that parseOutputForSystem expects
/// (production OGS always provides at least one reactant).
PhreeqcIOTestDriver makeNaClHaliteDriver(int num_systems, int num_threads = 1)
{
    PhreeqcIOTestDriver driver(num_systems, num_threads);
    driver.addComponent("Na");
    driver.addComponent("Cl");
    driver.setChargeBalance(ChargeBalance::Unspecified);
    driver.addEquilibriumPhase("Halite", /*saturation_index=*/-50.0,
                               /*initial_molality=*/0.0);
    driver.build();
    return driver;
}

/// Simple Calcite dissolution rate.  Note that `ReactionRate::operator<<`
/// already prepends a BASIC line number, so the statements must NOT carry
/// leading line numbers.  Avoids GOTO because the renumbered labels would not
/// match the original positions.
std::vector<std::string> calciteDissolutionRate()
{
    return {"si_cc = SI(\"Calcite\")", "rate = 1e-6 * (1 - 10^si_cc)",
            "IF (rate < 0) THEN rate = 0", "moles = rate * TIME", "SAVE moles"};
}
}  // namespace

// ============================================================================
// Group A: Dump Serialization
// ============================================================================

TEST(PhreeqcDump, PrintProducesDumpDirective)
{
    Dump dump("test_phreeqc_dump");
    std::ostringstream os;
    dump.print(os, 5);

    std::string const result = os.str();
    EXPECT_THAT(result, testing::HasSubstr("DUMP"));
    EXPECT_THAT(result, testing::HasSubstr("-solution 1-5"));
    EXPECT_THAT(result, testing::HasSubstr("-file"));
    EXPECT_THAT(result, testing::HasSubstr("END"));
}

TEST(PhreeqcDump, ReadDumpFromStringParsesMultipleSolutions)
{
    // Two SOLUTION_RAW blocks, each ending with -gammas
    std::string const dump_content =
        "SOLUTION_RAW 1\n"
        "  -temp 25\n"
        "  -gammas\n"
        "SOLUTION_RAW 2\n"
        "  -temp 30\n"
        "  -gammas\n";

    Dump dump;

    // num_systems refers to the total number of reactive transport
    // cells in the simulation.
    // PHREEQC uses IDs 1..N for the current timestep's solutions and
    // (N+1)..2N for the previous timestep's state loaded from the dump.
    // num_systems is intentionally different from the 2 solutions in
    // dump_content to verify that ID remapping uses num_systems,
    // not the solution count.
    std::size_t const num_systems = 3;
    dump.readDumpFromString(dump_content, num_systems);

    ASSERT_EQ(dump.aqueous_solutions_prev.size(), 2u);
    // IDs should be remapped: num_systems + 0 + 1 = 4, num_systems + 1 + 1 = 5
    EXPECT_THAT(dump.aqueous_solutions_prev[0],
                testing::HasSubstr("SOLUTION_RAW 4"));
    EXPECT_THAT(dump.aqueous_solutions_prev[1],
                testing::HasSubstr("SOLUTION_RAW 5"));
}

TEST(PhreeqcDump, ReadDumpFromStringForSystemParallelMode)
{
    // In parallel mode, each PHREEQC instance solves a single cell and
    // produces a dump containing SOLUTION_RAW 1. readDumpFromStringForSystem
    // places that result into the correct slot (system_id) of the
    // aqueous_solutions_prev vector and remaps the ID to
    // num_systems + system_id + 1 to avoid colliding with the current
    // timestep's solutions (IDs 1..N).
    std::string const dump_content =
        "SOLUTION_RAW 1\n"
        "  -temp 25\n"
        "  -gammas\n";

    Dump dump;
    std::size_t const num_systems = 4;
    // pick an arbitrary choice to keep an index that's not 0 or the last slot.
    std::size_t const system_id = 2;
    dump.readDumpFromStringForSystem(dump_content, system_id, num_systems);

    ASSERT_EQ(dump.aqueous_solutions_prev.size(), num_systems);
    // ID remapped: num_systems + system_id + 1 = 7
    EXPECT_THAT(dump.aqueous_solutions_prev[system_id],
                testing::HasSubstr("SOLUTION_RAW 7"));
    // Other slots should be empty
    EXPECT_TRUE(dump.aqueous_solutions_prev[0].empty());
    EXPECT_TRUE(dump.aqueous_solutions_prev[1].empty());
    EXPECT_TRUE(dump.aqueous_solutions_prev[3].empty());
}

TEST(PhreeqcDump, ReadDumpFromStringForSystemMultipleSystems)
{
    Dump dump;
    std::size_t const num_systems = 3;

    for (std::size_t i = 0; i < num_systems; ++i)
    {
        std::string const content =
            "SOLUTION_RAW 1\n"
            "  -temp " +
            std::to_string(25.0 + static_cast<double>(i)) +
            "\n"
            "  -gammas\n";
        dump.readDumpFromStringForSystem(content, i, num_systems);
    }

    ASSERT_EQ(dump.aqueous_solutions_prev.size(), num_systems);
    for (std::size_t i = 0; i < num_systems; ++i)
    {
        std::string const expected_id =
            "SOLUTION_RAW " + std::to_string(num_systems + i + 1);
        EXPECT_THAT(dump.aqueous_solutions_prev[i],
                    testing::HasSubstr(expected_id))
            << "System " << i << " should have ID " << expected_id;
        EXPECT_FALSE(dump.aqueous_solutions_prev[i].empty())
            << "System " << i << " should not be empty";
    }
}

TEST(PhreeqcDump, ReadDumpFileParsesStreamInput)
{
    std::string const dump_content =
        "SOLUTION_RAW 1\n"
        "  -temp 25\n"
        "  -gammas\n"
        "SOLUTION_RAW 2\n"
        "  -temp 30\n"
        "  -gammas\n";

    std::size_t const num_systems = 3;

    // Parse via readDumpFromString
    Dump dump_string;
    dump_string.readDumpFromString(dump_content, num_systems);

    // Parse via readDumpFile (istream overload)
    Dump dump_file;
    std::istringstream iss(dump_content);
    dump_file.readDumpFile(iss, num_systems);

    ASSERT_EQ(dump_string.aqueous_solutions_prev.size(),
              dump_file.aqueous_solutions_prev.size());
    for (std::size_t i = 0; i < dump_string.aqueous_solutions_prev.size(); ++i)
    {
        EXPECT_EQ(dump_string.aqueous_solutions_prev[i],
                  dump_file.aqueous_solutions_prev[i])
            << "Mismatch at solution " << i;
    }
}

TEST(PhreeqcDump, ReadDumpFromStringHandlesCarriageReturn)
{
    // Same content but with \r\n line endings
    std::string const dump_content =
        "SOLUTION_RAW 1\r\n"
        "  -temp 25\r\n"
        "  -gammas\r\n";

    Dump dump_cr;
    dump_cr.readDumpFromString(dump_content, 2);

    Dump dump_lf;
    dump_lf.readDumpFromString(
        "SOLUTION_RAW 1\n"
        "  -temp 25\n"
        "  -gammas\n",
        2);

    ASSERT_EQ(dump_cr.aqueous_solutions_prev.size(), 1u);
    ASSERT_EQ(dump_lf.aqueous_solutions_prev.size(), 1u);
    EXPECT_EQ(dump_cr.aqueous_solutions_prev[0],
              dump_lf.aqueous_solutions_prev[0]);
}

// ============================================================================
// Group B: PhreeqcInstancePool
// ============================================================================

TEST(PhreeqcInstancePool, CreateSingleInstance)
{
    PhreeqcInstancePool pool(phreeqc_database, 1);
    EXPECT_EQ(pool.size(), 1);
    EXPECT_GE(pool.getInstanceForThread(0), 0);
}

TEST(PhreeqcInstancePool, CreateMultipleInstances)
{
    PhreeqcInstancePool pool(phreeqc_database, 4);
    EXPECT_EQ(pool.size(), 4);

    // All IDs should be distinct
    std::vector<int> ids;
    for (int i = 0; i < 4; ++i)
    {
        ids.push_back(pool.getInstanceForThread(i));
    }
    std::sort(ids.begin(), ids.end());
    auto last = std::unique(ids.begin(), ids.end());
    EXPECT_EQ(std::distance(ids.begin(), last), 4)
        << "All instance IDs should be distinct";
}

// Regression: PhreeqcIO must use the IPhreeqc id returned by createInstance,
// not a hard-coded 0.  Keep a foreign IPhreeqc alive while PhreeqcIO is
// constructed and exercised; the foreign instance must remain operational
// before, during, and after PhreeqcIO's lifetime.
TEST(PhreeqcInstancePool, PhreeqcIOWorksWithForeignInstanceAlive)
{
    // The pool owns destruction of foreign_id, so no separate RAII guard is
    // needed.
    PhreeqcInstancePool foreign_pool(phreeqc_database, /*num_instances=*/1);
    int const foreign_id = foreign_pool.getInstanceForThread(0);
    ASSERT_GE(foreign_id, 0);
    PhreeqcInstancePool::configureForStringIO(foreign_id);
    ASSERT_EQ(RunString(foreign_id, "SOLUTION 1\nEND\n"), IPQ_OK);

    {
        auto driver =
            makeNaClHaliteDriver(/*num_systems=*/1, /*num_threads=*/1);

        driver.setInitialState(0, {{"Na", 0.001}, {"Cl", 0.001}}, /*pH=*/7.0);
        EXPECT_NO_THROW(driver.run(/*dt=*/1.0));

        // Round-trip check: dilute NaCl at pH 7 is already in equilibrium,
        // so PHREEQC must leave pH untouched.  Trivial chemistry keeps a
        // failure pointing at the plumbing rather than at the solver.
        EXPECT_NEAR(driver.pH(0), 7.0, 0.05);
    }

    // PhreeqcIO destruction must not have killed the foreign instance.
    EXPECT_EQ(RunString(foreign_id, "SOLUTION 1\nEND\n"), IPQ_OK);
}

// ============================================================================
// Group C: OpenMP Threading
// ============================================================================

TEST(PhreeqcParallel, MultiThreadExecution)
{
    // Trivial NaCl chemistry -- pH stays at 7 by design so that a failure
    // points at the parallel-pool plumbing rather than at PHREEQC.  Each
    // system is set up with a different Na/Cl total to ensure threads do
    // not silently share state.
#ifdef _OPENMP
    // Do not call omp_set_num_threads here: that would persist across tests
    // and silently change later tests' parallelism.  Inherit whatever the
    // environment provides and skip if it is not enough to actually exercise
    // the multi-thread path.
    int const num_threads = std::max(omp_get_max_threads(), 1);
    if (num_threads < 2)
    {
        GTEST_SKIP() << "MultiThreadExecution needs >= 2 OpenMP threads "
                        "(omp_get_max_threads() = "
                     << num_threads << ").";
    }
#else
    GTEST_SKIP() << "Built without OpenMP; multi-thread plumbing not under "
                    "test.";
    int const num_threads = 1;
#endif

    int const num_systems = 8;
    auto driver = makeNaClHaliteDriver(num_systems, num_threads);

    for (int i = 0; i < num_systems; ++i)
    {
        double const na = 0.001 * (i + 1);
        driver.setInitialState(i, {{"Na", na}, {"Cl", na}}, /*pH=*/7.0);
    }

    // One call drives the OpenMP loop inside
    // PhreeqcIO::executeSpeciationCalculationParallel for all systems.
    driver.run(/*dt=*/1.0);

    for (int i = 0; i < num_systems; ++i)
    {
        EXPECT_NEAR(driver.pH(i), 7.0, 0.05) << "System " << i;
    }
}

// ============================================================================
// Group D: PHREEQC Round-Trip
// ============================================================================

TEST(PhreeqcRoundTrip, ExchangeRoundTrip)
{
    // PhreeqcIO with an ExchangeSite triggers Dump allocation in the driver,
    // so writeSystemBlock emits a SOLUTION_RAW + EXCHANGE -equilibrate-with-
    // prev block from the second run onward.  We run twice to exercise that
    // dump path and check the chemistry stays bounded.
    PhreeqcIOTestDriver driver(/*num_systems=*/1, /*num_threads=*/1);
    driver.addComponent("Na").addComponent("Ca").addComponent("Cl");
    driver.setChargeBalance(ChargeBalance::pH);
    driver.addExchanger("X", /*initial_molality=*/1.1e-3);
    driver.build();

    driver.setInitialState(0, {{"Na", 0.001}, {"Ca", 0.0005}, {"Cl", 0.002}},
                           /*pH=*/7.0);

    driver.run(/*dt=*/1.0);

    // PHREEQC's charge-balanced pH for this Na/Ca/Cl solution equilibrated with
    // the X exchanger (phreeqc.dat).  The exchanger adsorbs Ca/Na and the pH
    // settles just below neutral.
    double const ph_step1 = driver.pH(0);
    EXPECT_NEAR(ph_step1, 6.995, 1e-3);

    // Second call uses the dump captured by run 1 (Dump
    // ::aqueous_solutions_prev was populated inside
    // executeSpeciationCalculationParallel).  With unchanged inputs it must
    // reproduce step 1 exactly; a loose band would hide dump-readback drift.
    driver.run(/*dt=*/1.0);
    double const ph_step2 = driver.pH(0);
    EXPECT_NEAR(ph_step2, ph_step1, 1e-9);
}

TEST(PhreeqcRoundTrip, EquilibriumPhasesRoundTrip)
{
    // Calcite at SI=0 buffers pH above neutral.  Single PhreeqcIO call drives
    // writeSystemBlock to emit an EQUILIBRIUM_PHASES block and
    // parseOutputForSystem to write the buffered pH back.
    PhreeqcIOTestDriver driver(/*num_systems=*/1, /*num_threads=*/1);
    driver.addComponent("Ca").addComponent("C(4)", "HCO3").addComponent("Cl");
    driver.setChargeBalance(ChargeBalance::pH);
    driver.addEquilibriumPhase("Calcite", /*saturation_index=*/0.0,
                               /*initial_molality=*/0.01);
    driver.build();

    driver.setInitialState(0, {{"Ca", 0.001}, {"C(4)", 0.002}, {"Cl", 0.001}},
                           /*pH=*/7.0);

    driver.run(/*dt=*/1.0);

    // Calcite-buffered pH computed by PHREEQC for this Ca/C(4)/Cl solution
    // (phreeqc.dat); equilibration with Calcite at SI=0 raises it above
    // neutral.
    double const ph_after = driver.pH(0);
    EXPECT_GT(ph_after, 7.0) << "Calcite-buffered water should have pH > 7";
    EXPECT_NEAR(ph_after, 7.372, 1e-3);
}

TEST(PhreeqcRoundTrip, KineticsRoundTrip)
{
    // Calcite dissolution kinetics over one hour.  The driver passes the
    // matching ReactionRate script into PhreeqcIO's constructor; PhreeqcIO
    // emits both a RATES header (writeInputHeader) and a KINETICS block
    // (writeSystemBlock).  We verify Ca rises and pH increases as Calcite
    // dissolves.
    PhreeqcIOTestDriver driver(/*num_systems=*/1, /*num_threads=*/1);
    driver.addComponent("Ca").addComponent("C(4)", "HCO3").addComponent("Cl");
    driver.setChargeBalance(ChargeBalance::pH);
    driver.addKineticPhase("Calcite", /*chemical_formula=*/"CaCO3",
                           /*initial_molality=*/0.01, calciteDissolutionRate());
    driver.build();

    driver.setInitialState(0, {{"Ca", 0.0005}, {"C(4)", 0.001}, {"Cl", 0.001}},
                           /*pH=*/7.0);

    driver.run(/*dt=*/3600.0);

    // PHREEQC-computed pH after one hour of Calcite dissolution kinetics
    // (phreeqc.dat); dissolving Calcite consumes H+ and raises pH above
    // neutral.
    double const ph_after = driver.pH(0);
    EXPECT_GT(ph_after, 7.0) << "Calcite dissolution should raise pH";
    EXPECT_NEAR(ph_after, 7.603, 1e-3);

    // Ca total inventory rises from its starting value of 0.5 mmol/kgw to
    // PHREEQC-computed 1.416 mmol/kgw as Calcite dissolves over one hour
    // (phreeqc.dat, kinetic rate calciteDissolutionRate()).
    double const ca_after = driver.getConcentration(/*component_id=*/0,
                                                    /*chemical_system_id=*/0);
    EXPECT_GE(ca_after, 0.0005);
    EXPECT_NEAR(ca_after, 1.4164e-3, 1e-6);
}

TEST(PhreeqcRoundTrip, CombinedReactionsRoundTrip)
{
    PhreeqcIOTestDriver driver(/*num_systems=*/1, /*num_threads=*/1);
    driver.addComponent("Na")
        .addComponent("Ca")
        .addComponent("C(4)", "HCO3")
        .addComponent("Cl");
    driver.setChargeBalance(ChargeBalance::pH);
    driver.addEquilibriumPhase("Calcite", 0.0, 0.01);
    driver.addExchanger("X", 1.1e-3);
    driver.build();

    driver.setInitialState(
        0,
        {{"Na", 1.0e-3}, {"Ca", 0.5e-3}, {"C(4)", 2.0e-3}, {"Cl", 2.0e-3}},
        /*pH=*/7.0);

    driver.run(/*dt=*/1.0);

    // PHREEQC-computed pH for the combined Calcite + X-exchange system
    // (phreeqc.dat).
    double const ph_step1 = driver.pH(0);
    EXPECT_NEAR(ph_step1, 7.225, 1e-3);

    // Second call exercises writeSystemBlock's dump-restore branch on top of
    // equilibrium phases and exchange.  With unchanged inputs it must reproduce
    // step 1 exactly; a loose band would hide dump-readback drift.
    driver.run(/*dt=*/1.0);
    double const ph_step2 = driver.pH(0);
    EXPECT_NEAR(ph_step2, ph_step1, 1e-9);
}

TEST(PhreeqcRoundTrip, ParallelExchangeRoundTrip)
{
    int const num_systems = 4;

#ifdef _OPENMP
    int const num_threads = std::max(omp_get_max_threads(), 1);
#else
    int const num_threads = 1;
#endif

    PhreeqcIOTestDriver driver(num_systems, num_threads);
    driver.addComponent("Na").addComponent("Ca").addComponent("Cl");
    driver.setChargeBalance(ChargeBalance::pH);
    driver.addExchanger("X", 1.1e-3);
    driver.build();

    for (int i = 0; i < num_systems; ++i)
    {
        double const na = 0.001 * (i + 1);
        double const ca = 0.0005 * (i + 1);
        double const cl = na + 2 * ca;
        driver.setInitialState(i, {{"Na", na}, {"Ca", ca}, {"Cl", cl}},
                               /*pH=*/7.0);
    }

    // First call: dump is empty; writeSystemBlock equilibrates exchange with
    // the current solution.  Second call: dump from step 1 drives the
    // -equilibrate-with-prev branch and the parallel dump-readback path in
    // executeSpeciationCalculationParallel.
    // PHREEQC-computed step-1 pH per system (phreeqc.dat); pH drops slightly as
    // the Na/Ca load rises with system index.
    std::array<double, 4> const expected_pH_step1 = {6.99501, 6.99301, 6.99120,
                                                     6.98952};

    driver.run(/*dt=*/1.0);
    std::vector<double> pH_step1(num_systems);
    for (int i = 0; i < num_systems; ++i)
    {
        pH_step1[i] = driver.pH(i);
        EXPECT_NEAR(pH_step1[i], expected_pH_step1[i], 1e-3)
            << "System " << i << " step 1";
    }

    // Without changing inputs, step 2 must reproduce step 1 closely.  This
    // catches dump-readback bugs that would only show up as drift, not as a
    // value outside the loose 4..10 band.
    driver.run(/*dt=*/1.0);
    for (int i = 0; i < num_systems; ++i)
    {
        EXPECT_NEAR(driver.pH(i), pH_step1[i], 1e-6)
            << "System " << i << ": step-2 pH drifted from step-1";
    }
}

// ============================================================================
// Group E: MPI Error Propagation
// ============================================================================

#ifdef USE_PETSC
struct MPI_PhreeqcParallel : public ::testing::Test
{
    BaseLib::MPI::Mpi mpi{};
};

TEST_F(MPI_PhreeqcParallel, AllRanksSucceed)
{
    // Each rank runs an independent PhreeqcIO over a one-cell mesh with
    // slightly different concentrations.  PhreeqcIO::
    // executeSpeciationCalculationParallel already calls allRanksThrowOrNone
    // internally (PhreeqcIO.cpp:1277); this test simply verifies the call
    // returns without throwing on all ranks.
    auto driver = makeNaClHaliteDriver(/*num_systems=*/1, /*num_threads=*/1);

    double const c = 0.001 * (mpi.rank + 1);
    driver.setInitialState(0, {{"Na", c}, {"Cl", c}}, /*pH=*/7.0);

    EXPECT_NO_THROW(driver.run(/*dt=*/1.0));
    EXPECT_NEAR(driver.pH(0), 7.0, 0.05);
}

TEST_F(MPI_PhreeqcParallel, ErrorPropagation)
{
    if (mpi.size < 2)
    {
        GTEST_SKIP() << "Requires at least 2 MPI ranks.";
    }

    // Rank 1 declares a nonexistent mineral as an EquilibriumReactant.
    // When PhreeqcIO::writeSystemBlock emits the EQUILIBRIUM_PHASES block,
    // PHREEQC will reject the unknown mineral name, RunString returns
    // failure, executeSpeciationCalculationParallel collects the failure in
    // failed_systems, and the internal allRanksThrowOrNone propagates the
    // error to every rank.  We confirm run() throws on every rank.
    PhreeqcIOTestDriver driver(/*num_systems=*/1, /*num_threads=*/1);
    driver.addComponent("Na").addComponent("Cl");
    driver.setChargeBalance(ChargeBalance::Unspecified);
    // Every rank needs at least one valid Solid constituent, otherwise the
    // test Medium would carry an empty <Solid> phase and CreatePhase would
    // OGS_FATAL during build() before any chemistry runs.  Halite at SI=-50
    // with mol=0 never precipitates, so it is an inert placeholder (see
    // makeNaClHaliteDriver).
    driver.addEquilibriumPhase("Halite", /*saturation_index=*/-50.0,
                               /*initial_molality=*/0.0);
    if (mpi.rank == 1)
    {
        driver.addEquilibriumPhase("NonexistentMineral",
                                   /*saturation_index=*/0.0,
                                   /*initial_molality=*/1.0);
    }
    driver.build();

    driver.setInitialState(0, {{"Na", 0.001}, {"Cl", 0.001}}, /*pH=*/7.0);

    EXPECT_THROW(driver.run(/*dt=*/1.0), std::exception);
}
#endif
