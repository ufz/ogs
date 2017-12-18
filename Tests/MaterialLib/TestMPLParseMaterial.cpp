/**
 * \author Norbert Grunwald
 * \date   Sep 22, 2017
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include <memory>
#include <sstream>

#include "Tests/TestTools.h"

#include "MaterialLib/MPL/mpMedium.h"

namespace MPL = MaterialPropertyLib;

//----------------------------------------------------------------------------
// Test density models.
MPL::Medium createTestMaterial(std::string const& xml)
{
    auto const ptree = readXml(xml);
    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& config =
        conf.getConfigSubtree("media").getConfigSubtree("medium");

    return MPL::Medium(config);
}
/// A simple structure for a string-based component used to create
/// an XML-Structure.
struct Component
{
    std::vector<std::string> property;
    Component() : property(MPL::number_of_property_enums) {}
};
/// A simple structure for a string-based phase used to create
/// an XML-Structure.
struct Phase
{
    std::vector<Component> component;
    std::vector<std::string> property;
    Phase(std::size_t componentNumber)
        : component(componentNumber), property(MPL::number_of_property_enums)
    {
    }
};
/// A simple structure for a string-based medium used to create
/// an XML-Structure
struct Medium
{
    std::vector<Phase> phases;
    std::vector<std::string> property;
    Medium(std::vector<std::size_t> topology)
        : property(MPL::number_of_property_enums)
    {
        for (auto p : topology)
            phases.push_back(Phase(p));
    }
};

/// Short method creating a number of blanks used for indentation.
std::string indent(std::size_t ind)
{
    return std::string(ind, ' ');
}

/// This method creates an XML-snippet of a sigle property.
std::string makeProperty(std::size_t ind, std::size_t index,
                         std::string property)
{
    std::stringstream sProperty;

    if (property == "")
        return "";

    sProperty << indent(ind) << "<property>\n";
    sProperty << indent(ind + 2) << "<name>" << MPL::convertEnumToString[index]
              << "</name>\n";

    // What follows is a C-style string-to-double-conversion...
    // seems so easy compared to c++-standard; Feel free to
    // change it to latest standard..
    if (double value = atof(property.c_str()))
    {
        sProperty << indent(ind + 2) << "<type>constant</type>\n";
        sProperty << indent(ind + 2) << "<value>" << value << "</value>\n";
    }
    else
        sProperty << indent(ind + 2) << "<type>" << property << "</type>\n";

    sProperty << indent(ind) << "</property>\n";
    return sProperty.str();
}

/// This method creates an XML-snippet of a property section.
std::string makeProperties(std::size_t ind, std::vector<std::string> properties)
{
    std::stringstream sProperties;
    std::vector<std::string> temp;
    bool empty(true);

    for (std::size_t p = 0; p < properties.size(); ++p)
    {
        std::string property = makeProperty(ind + 2, p, properties[p]);
        if ((property != "") && (p != static_cast<std::size_t>(MPL::name)))
        {
            temp.push_back(property);
            empty = false;
        }
        // If a property is not given, a default takes its place. Note
        // also that 'name' is an exceptional property, since it determines
        // specific components (thus it is defined outside of the properties-
        // tag).
    }

    if (empty)
        return "";
    sProperties << indent(ind) << "<properties>\n";
    for (auto str : temp)
        sProperties << str;
    sProperties << indent(ind) << "</properties>\n";
    return sProperties.str();
}

/// This method creates an XML-snippet of a single component.
std::string makeComponent(Component c)
{
    std::string componentProperties = makeProperties(14, c.property);
    std::stringstream component;

    component << indent(12) << "<component>\n";
    if (c.property[MPL::name] != "")
        component << indent(14) << "<name>" << c.property[MPL::name]
                  << "</name>\n";
    component << componentProperties;
    component << indent(12) << "</component>\n";
    return component.str();
}

/// This method creates an XML-snippet of the phase
/// components.
std::string makeComponents(std::vector<Component> components)
{
    std::stringstream sComponents;
    sComponents << indent(10) << "<components>\n";
    for (auto c : components)
    {
        std::string component = makeComponent(c);
        sComponents << component;
    }
    sComponents << indent(10) << "</components>\n";
    return sComponents.str();
}

/// This method creates an XML-snippet of a sigle material
/// phase.
std::string makePhase(Phase p)
{
    std::string phaseComponents = makeComponents(p.component);
    std::string phaseProperties = makeProperties(8, p.property);
    std::stringstream phase;

    phase << indent(8) << "<phase>\n";
    if (p.property[MPL::name] != "")
        phase << indent(10) << "<name>" << p.property[MPL::name] << "</name>\n";
    phase << phaseComponents;
    phase << phaseProperties;
    phase << indent(8) << "</phase>\n";
    return phase.str();
}

/// This method creates an XML-snippet of the material
/// phases.
std::string makePhases(std::vector<Phase> phases)
{
    std::stringstream sPhases;

    sPhases << indent(6) << "<phases>\n";
    for (auto p : phases)
    {
        std::string phase = makePhase(p);
        sPhases << phase;
    }
    sPhases << indent(6) << "</phases>\n";
    return sPhases.str();
}

/// This method creates the entire XML-tree structure from the string-
/// based medium specification object. I know, indentation is not
/// necessary for this, but my OCD kicked in...
std::string makeMedium(Medium m)
{
    std::string mediumPhases = makePhases(m.phases);
    std::string mediumProperties = makeProperties(6, m.property);
    std::stringstream medium;
    medium << "<media>\n";
    medium << indent(2) << "<medium>\n";
    if (m.property[MPL::name] != "")
        medium << indent(4) << "<name>" << m.property[MPL::name] << "</name>\n";
    medium << mediumPhases;
    medium << mediumProperties;
    medium << indent(2) << "</medium>\n";
    medium << "</media>\n";
    return medium.str();
}

/// A method used to obtain the name of a medium, phase, or component of a
/// material or of a specifier and to store them in two vectors for later
/// comparison.
void getNames(MPL::Property const& observation, std::string expectation,
              std::string defaultName, std::vector<std::string>* obs,
              std::vector<std::string>* exp)
{
    obs->push_back(MPL::getString(observation));
    if (expectation == "")
        exp->push_back(defaultName);
    else
        exp->push_back(expectation);
}

// This test represents an invariant test. First, several phases,
// components, and properties are generated (more or less randomly).
// Then, an XML-tree is generated from those information. The XML
// tree is then parsed into a configTree, from which an
// MPL::MaterialObject is generated.
// Last step compares the names (as well as the topology) of the
// Material object with the specified parameters.
TEST(Material, parseMaterials)
{
    // This is the topology of our new material: The size of the
    // topology vector determines the number of phases, while each
    // vector component refers to the number of components of that
    // phase.
    // The number of properties is fixed in each case and is
    // determined by the size of the PropertyEnum enumerator.
    std::vector<std::size_t> const mediumTopology = {1, 3, 2};

    // make a string from every property enumerator
    std::vector<std::string> property = MPL::convertEnumToString;

    Medium medium(mediumTopology);

    // the omnivagant medium:
    medium.property[MPL::name] = "luminiferous_aether";

    medium.phases[0].property[MPL::name] = "solid";
    medium.phases[1].property[MPL::name] = "liquid";
    medium.phases[2].property[MPL::name] = "gas";

    medium.phases[0].component[0].property[MPL::density] = "LinearTemperature";
    medium.phases[0].component[0].property[MPL::thermal_conductivity] = "0.654";
    medium.phases[0].component[0].property[MPL::reference_temperature] = "333";
    medium.phases[0].component[0].property[MPL::reference_density] = "2100.0";
    medium.phases[0].component[0].property[MPL::drhodT] = "-0.4";
    medium.phases[0].component[0].property[MPL::effective_stress] =
        "LinearEpsilon";

    medium.phases[1].component[0].property[MPL::name] = "Water";
    medium.phases[1].component[1].property[MPL::name] = "CarbonDioxide";
    medium.phases[1].component[2].property[MPL::name] = "SodiumChloride";

    medium.phases[1].property[MPL::density] = "Duan_2012";
    medium.phases[1].property[MPL::viscosity] = "Islam_Carlson_2012";

    medium.phases[2].component[0].property[MPL::name] = "Water";
    medium.phases[2].component[0].property[MPL::viscosity] = "IAPWS_2008";
    medium.phases[2].component[1].property[MPL::name] = "CarbonDioxide";
    medium.phases[2].component[1].property[MPL::viscosity] = "Fenghour_1998";

    medium.phases[2].property[MPL::density] = "Peng_Robinson_1976";
    medium.phases[2].property[MPL::viscosity] = "Buddenberg_Wilke_1949";

    medium.property[MPL::thermal_conductivity] = "AverageVolumeFraction";
    medium.property[MPL::saturation] = "Brooks_Corey_1964";
    medium.property[MPL::permeability] = "1.0e-12";
    medium.property[MPL::heat_capacity] = "AverageVolumeFraction";
    medium.property[MPL::relative_permeability] = "Mualem_1978";

    // create an actual MaterialProperty-Medium out of the specifier object
    auto m = createTestMaterial(makeMedium(medium));

    // those two vectors will actually be compared
    std::vector<std::string> expected;
    std::vector<std::string> observed;

    // get the names of the specifier and that of the medium and store
    // them in two vectors. If a name of the medium, one of the phases
    // or components is not specified, it is automatically replaced by
    // the default (which has to be the same as the default of the MPL
    getNames(m.property(MPL::name), medium.property[MPL::name], "no_name",
             &observed, &expected);

    // now we roam through all phases and components, finding their names
    // and storing them in the two vectors
    for (std::size_t p = 0; p < m.numberOfPhases(); ++p)
    {
        const auto& phase = m.phase(p);
        getNames(phase.property(MPL::name),
                 medium.phases[p].property[MPL::name], "no_name", &observed,
                 &expected);

        for (std::size_t c = 0; c < phase.numberOfComponents(); ++c)
        {
            const auto& component = phase.component(c);
            getNames(component.property(MPL::name),
                     medium.phases[p].component[c].property[MPL::name],
                     "no_name", &observed, &expected);
        }
    }

    // Now, the two vectors are compared. If there is some derivation,
    // we can easily locate the problem.
    ASSERT_EQ(expected, observed);
}
