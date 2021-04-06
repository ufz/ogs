/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 22, 2017
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include <sstream>

#include "MaterialLib/MPL/Medium.h"
#include "TestMPL.h"
#include "Tests/TestTools.h"

namespace MPL = MaterialPropertyLib;

/// A simple structure for a string-based component used to create
/// an XML-Structure.
struct Component
{
    std::vector<std::string> property;
    Component() : property(MPL::number_of_properties) {}
};
/// A simple structure for a string-based phase used to create
/// an XML-Structure.
struct Phase
{
    std::vector<Component> component;
    std::vector<std::string> property;
    explicit Phase(std::size_t const componentNumber)
        : component(componentNumber), property(MPL::number_of_properties)
    {
    }
};

/// A simple structure for a string-based medium used to create
/// an XML-Structure
struct Medium
{
    std::vector<Phase> phases;
    std::vector<std::string> property;
    explicit Medium(std::vector<std::size_t> const& phases_)
        : property(MPL::number_of_properties)
    {
        for (auto p : phases_)
        {
            phases.emplace_back(p);
        }
    }
};

/// Short method creating a number of blanks used for indentation.
std::string indent(std::size_t const index)
{
    return std::string(index, ' ');
}

/// This method creates an XML-snippet of a single property.
std::string makeProperty(std::size_t const ind, std::size_t const index,
                         std::string const& property)
{
    std::stringstream sProperty;

    if (property.empty())
    {
        return "";
    }

    sProperty << indent(ind) << "<property>\n";
    sProperty << indent(ind + 2) << "<name>"
              << MPL::property_enum_to_string[index] << "</name>\n";

    if (double const value = std::atof(property.c_str()))
    {
        sProperty << indent(ind + 2) << "<type>Constant</type>\n";
        sProperty << indent(ind + 2) << "<value>" << value << "</value>\n";
    }
    else
    {
        sProperty << indent(ind + 2) << "<type>" << property << "</type>\n";
    }

    sProperty << indent(ind) << "</property>\n";
    return sProperty.str();
}

/// This method creates an XML-snippet of a property section.
std::string makeProperties(std::size_t const ind,
                           std::vector<std::string> const& properties)
{
    std::vector<std::string> temp;

    for (std::size_t p = 0; p < properties.size(); ++p)
    {
        std::string property = makeProperty(ind + 2, p, properties[p]);
        if ((!property.empty()) && (p != static_cast<std::size_t>(MPL::name)))
        {
            temp.push_back(property);
        }
        // If a property is not given, a default takes its place. Note also that
        // 'name' is an exceptional property, since it determines specific
        // components (thus it is defined outside of the properties- tag).
    }

    if (temp.empty())
    {
        return "";
    }
    std::stringstream sProperties;
    sProperties << indent(ind) << "<properties>\n";
    for (auto const& str : temp)
    {
        sProperties << str;
    }
    sProperties << indent(ind) << "</properties>\n";
    return sProperties.str();
}

/// This method creates an XML-snippet of a single component.
std::string makeComponent(Component const& c)
{
    std::stringstream component;
    component << indent(12) << "<component>\n";
    component << indent(14) << "<name>" << c.property[MPL::name] << "</name>\n";
    component << makeProperties(14, c.property);
    component << indent(12) << "</component>\n";
    return component.str();
}

/// This method creates an XML-snippet of the phase components.
std::string makeComponents(std::vector<Component> const& components)
{
    std::stringstream sComponents;
    sComponents << indent(10) << "<components>\n";
    for (auto const& c : components)
    {
        sComponents << makeComponent(c);
    }
    sComponents << indent(10) << "</components>\n";
    return sComponents.str();
}

/// This method creates an XML-snippet of a single material phase.
std::string makePhase(Phase const& p)
{
    std::stringstream phase;

    phase << indent(8) << "<phase>\n";
    phase << indent(10) << "<type>" << p.property[MPL::name] << "</type>\n";
    phase << makeComponents(p.component);
    phase << makeProperties(8, p.property);
    phase << indent(8) << "</phase>\n";
    return phase.str();
}

/// This method creates an XML-snippet of the material phases.
std::string makePhases(std::vector<Phase> const& phases)
{
    std::stringstream sPhases;

    sPhases << indent(6) << "<phases>\n";
    for (auto const& p : phases)
    {
        sPhases << makePhase(p);
    }
    sPhases << indent(6) << "</phases>\n";
    return sPhases.str();
}

/// This method creates the entire XML-tree structure from the string- based
/// medium specification object. I know, indentation is not necessary for this,
/// but my OCD kicked in...
std::string makeMedium(Medium const& m)
{
    std::stringstream medium;
    medium << indent(2) << "<medium>\n";
    medium << makePhases(m.phases);
    medium << makeProperties(6, m.property);
    medium << indent(2) << "</medium>\n";
    return medium.str();
}

/// A method used to obtain the name of a medium, phase, or component of a
/// material or of a specifier and to store them in two vectors for later
/// comparison.
void getNames(std::string const& observed_name, std::string const& expectation,
              std::string const& defaultName, std::vector<std::string>& obs,
              std::vector<std::string>& exp)
{
    obs.push_back(observed_name);
    if (expectation.empty())
    {
        exp.push_back(defaultName);
    }
    else
    {
        exp.push_back(expectation);
    }
}

// This test represents an invariant test. First, several phases, components,
// and properties are generated (more or less randomly).
// Then, an XML-tree is generated from those information. The XML tree is then
// parsed into a configTree, from which an MPL::MaterialObject is generated.
// Last step compares the names (as well as the topology) of the Material object
// with the specified parameters.
TEST(Material, parseMaterials)
{
    // This is the topology of our new material: The size of the
    // topology vector determines the number of phases, while each
    // vector component refers to the number of components of that
    // phase.
    // The number of properties is fixed in each case and is
    // determined by the size of the PropertyEnum enumerator.
    std::vector<std::size_t> const mediumTopology = {1, 1, 1};

    Medium medium(mediumTopology);

    // the omnivagant medium:
    medium.property[MPL::name] = "luminiferous_aether";

    medium.phases[0].property[MPL::name] = "Solid";
    medium.phases[1].property[MPL::name] = "AqueousLiquid";
    medium.phases[2].property[MPL::name] = "Gas";

    medium.phases[0].component[0].property[MPL::thermal_conductivity] = "0.654";
    medium.phases[0].component[0].property[MPL::reference_temperature] = "333";
    medium.phases[0].component[0].property[MPL::reference_density] = "2100.0";
    medium.phases[0].component[0].property[MPL::drhodT] = "-0.4";

    medium.phases[0].component[0].property[MPL::name] = "VerySolid";
    medium.phases[1].component[0].property[MPL::name] = "Water";
    medium.phases[2].component[0].property[MPL::name] = "SuperFluid";
    medium.phases[2].component[0].property[MPL::thermal_conductivity] = "1";
    medium.property[MPL::permeability] = "1.0e-12";

    // create an actual MaterialProperty-Medium out of the specifier object
    auto const m = Tests::createTestMaterial(makeMedium(medium));

    // those two vectors will actually be compared
    std::vector<std::string> expected;
    std::vector<std::string> observed;

    // now we roam through all phases and components, finding their names
    // and storing them in the two vectors
    for (std::size_t p = 0; p < m->numberOfPhases(); ++p)
    {
        const auto& phase = m->phase(p);
        getNames(phase.name, medium.phases[p].property[MPL::name], "no_name",
                 observed, expected);

        for (std::size_t c = 0; c < phase.numberOfComponents(); ++c)
        {
            const auto& component = phase.component(c);
            getNames(component.name,
                     medium.phases[p].component[c].property[MPL::name],
                     "no_name", observed, expected);
        }
    }

    // Now, the two vectors are compared. If there is some derivation,
    // we can easily locate the problem.
    ASSERT_EQ(expected, observed);
}
