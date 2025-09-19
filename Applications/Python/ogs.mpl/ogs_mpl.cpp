/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <range/v3/range/conversion.hpp>
#include <range/v3/view/transform.hpp>

#include "BaseLib/Logging.h"
#include "BaseLib/MPI.h"
#include "MaterialLib/MPL/Properties/Constant.h"
#include "MaterialLib/MPL/Properties/Linear.h"
#include "MaterialLib/MPL/Property.h"
#include "MathLib/Point3d.h"
#include "ParameterLib/SpatialPosition.h"

namespace py = pybind11;

using namespace MaterialPropertyLib;
using namespace ParameterLib;

namespace
{
// Convert NumPy array to Point3d if not None. This function is there to avoid
// pulling the MeshLib::Point3d object into python, but use std::array (which
// has all automatic conversions in pybind11).
std::optional<MathLib::Point3d> pythonToSpatialPositionCoords(
    py::object const& coordinates)
{
    if (coordinates.is_none())
    {
        return std::nullopt;
    }
    return std::make_optional<MathLib::Point3d>(
        py::cast<std::array<double, 3>>(coordinates));
}

// Convert SpatialPositon to python array if coordinates are set. This function
// is there to avoid pulling the MeshLib::Point3d object into python, but use
// python arrays.
py::object spatialPositionCoordsToPython(SpatialPosition const& pos)
{
    auto const& opt_coords = pos.getCoordinates();
    if (opt_coords)
    {
        return py::array_t<double>(3, opt_coords->data());
    }
    return py::none();
}

void bindSpatialPosition(py::module_& m)
{
    py::class_<SpatialPosition>(m, "SpatialPosition",
                                R"pbdoc(
    Describes a spatial position within a mesh or domain.

    A SpatialPosition may refer to a node (via node ID), an element (via element ID),
    or a physical point in space (via coordinates). It is typically used when evaluating
    material properties that depend on the spatial context within a simulation.

    The position may be partially specified. For example:
    - Only coordinates for geometric evaluation
    - Only node ID for nodal properties
    - Only element ID for element-wise values
    - Or any combination of the above

    The class is used as a parameter to property evaluations in the OGS material properties library.
)pbdoc")
        .def(py::init(
                 [](std::optional<std::size_t> node_id,
                    std::optional<std::size_t>
                        element_id,
                    py::object const& coordinates)
                 {
                     return SpatialPosition(
                         node_id, element_id,
                         pythonToSpatialPositionCoords(coordinates));
                 }),
             py::arg("node_id") = std::nullopt,
             py::arg("element_id") = std::nullopt,
             py::arg("coords") = py::none(),
             R"pbdoc(
             SpatialPosition(node_id=None, element_id=None, coords=None)

             Parameters:
                 node_id (int, optional): Node ID
                 element_id (int, optional): Element ID
                 coords (array-like of 3 floats, optional): Coordinates
             )pbdoc")

        .def_property(
            "node_id",
            [](SpatialPosition const& pos) { return pos.getNodeID(); },
            [](SpatialPosition& pos, std::size_t const id)
            { pos.setNodeID(id); },
            R"pbdoc(
            Node ID of the spatial position.

            This property can be read and set from Python. Setting to None is not supported.
            )pbdoc")

        .def_property(
            "element_id",
            [](SpatialPosition const& pos) { return pos.getElementID(); },
            [](SpatialPosition& pos, std::size_t const id)
            { pos.setElementID(id); },
            R"pbdoc(
            Element ID of the spatial position.

            This property can be read and set from Python. Setting to None is not supported.
            )pbdoc")

        .def_property(
            "coordinates",
            spatialPositionCoordsToPython,
            [](SpatialPosition& pos, py::object const& coordinates)
            {
                pos.setCoordinates(
                    pythonToSpatialPositionCoords(coordinates).value());
            },
            R"pbdoc(
            Coordinates of the spatial position as a 3-element array.

            This property can be read and set from Python. Use a 3-element list, tuple,
            or NumPy array. Setting to None is not supported.
            )pbdoc")

        .def(
            "__repr__",
            [](SpatialPosition const& pos)
            {
                auto const node_id =
                    pos.getNodeID() ? std::to_string(*pos.getNodeID()) : "None";
                auto const element_id =
                    pos.getElementID() ? std::to_string(*pos.getElementID())
                                       : "None";
                auto const coords = spatialPositionCoordsToPython(pos);

                return "<SpatialPosition(" + node_id + ", " + element_id +
                       ", " + py::str(coords).cast<std::string>() + ")>";
            },
            R"pbdoc(
            Return string representation of the SpatialPosition.
            )pbdoc");
}

void bindVariableArray(py::module_& m)
{
    py::class_<VariableArray>(m, "VariableArray")
        .def(py::init<>())
        .def_readwrite("capillary_pressure", &VariableArray::capillary_pressure)
        .def_readwrite("concentration", &VariableArray::concentration)
        .def_readwrite("deformation_gradient",
                       &VariableArray::deformation_gradient)
        .def_readwrite("density", &VariableArray::density)
        .def_readwrite("effective_pore_pressure",
                       &VariableArray::effective_pore_pressure)
        .def_readwrite("enthalpy", &VariableArray::enthalpy)
        .def_readwrite("enthalpy_of_evaporation",
                       &VariableArray::enthalpy_of_evaporation)
        .def_readwrite("equivalent_plastic_strain",
                       &VariableArray::equivalent_plastic_strain)
        .def_readwrite("fracture_aperture", &VariableArray::fracture_aperture)
        .def_readwrite("grain_compressibility",
                       &VariableArray::grain_compressibility)
        .def_readwrite("liquid_phase_pressure",
                       &VariableArray::liquid_phase_pressure)
        .def_readwrite("liquid_saturation", &VariableArray::liquid_saturation)
        .def_readwrite("mechanical_strain", &VariableArray::mechanical_strain)
        .def_readwrite("molar_mass", &VariableArray::molar_mass)
        .def_readwrite("molar_mass_derivative",
                       &VariableArray::molar_mass_derivative)
        .def_readwrite("molar_fraction", &VariableArray::molar_fraction)
        .def_readwrite("gas_phase_pressure", &VariableArray::gas_phase_pressure)
        .def_readwrite("porosity", &VariableArray::porosity)
        .def_readwrite("solid_grain_pressure",
                       &VariableArray::solid_grain_pressure)
        .def_readwrite("stress", &VariableArray::stress)
        .def_readwrite("temperature", &VariableArray::temperature)
        .def_readwrite("total_strain", &VariableArray::total_strain)
        .def_readwrite("total_stress", &VariableArray::total_stress)
        .def_readwrite("transport_porosity", &VariableArray::transport_porosity)
        .def_readwrite("vapour_pressure", &VariableArray::vapour_pressure)
        .def_readwrite("volumetric_strain", &VariableArray::volumetric_strain);
}

void bindVariableEnum(py::module_& m)
{
    py::enum_<Variable>(m, "Variable")
        .value("capillary_pressure", Variable::capillary_pressure)
        .value("concentration", Variable::concentration)
        .value("deformation_gradient", Variable::deformation_gradient)
        .value("density", Variable::density)
        .value("effective_pore_pressure", Variable::effective_pore_pressure)
        .value("enthalpy", Variable::enthalpy)
        .value("enthalpy_of_evaporation", Variable::enthalpy_of_evaporation)
        .value("equivalent_plastic_strain", Variable::equivalent_plastic_strain)
        .value("fracture_aperture", Variable::fracture_aperture)
        .value("grain_compressibility", Variable::grain_compressibility)
        .value("liquid_phase_pressure", Variable::liquid_phase_pressure)
        .value("liquid_saturation", Variable::liquid_saturation)
        .value("mechanical_strain", Variable::mechanical_strain)
        .value("molar_mass", Variable::molar_mass)
        .value("molar_mass_derivative", Variable::molar_mass_derivative)
        .value("molar_fraction", Variable::molar_fraction)
        .value("gas_phase_pressure", Variable::gas_phase_pressure)
        .value("porosity", Variable::porosity)
        .value("solid_grain_pressure", Variable::solid_grain_pressure)
        .value("stress", Variable::stress)
        .value("temperature", Variable::temperature)
        .value("total_strain", Variable::total_strain)
        .value("total_stress", Variable::total_stress)
        .value("transport_porosity", Variable::transport_porosity)
        .value("vapour_pressure", Variable::vapour_pressure)
        .value("volumetric_strain", Variable::volumetric_strain)
        .export_values();
}

void bindProperty(py::module_& m)
{
    py::class_<Property>(m, "Property", R"pbdoc(
        Base class for material properties.
    )pbdoc")
        .def(
            "value",
            [](const Property& p, const VariableArray& va,
               const SpatialPosition& pos, double t, double dt)
            { return p.value(va, pos, t, dt); },
            py::arg("variable_array"),
            py::arg("pos"),
            py::arg("t"),
            py::arg("dt"),
            R"pbdoc(
             Evaluate the property value.

             Parameters:
                 variable_array: Current time step values.
                 pos: Spatial position.
                 t: Current time.
                 dt: Time step size.
             )pbdoc")

        .def(
            "value",
            [](const Property& p, const VariableArray& va,
               const VariableArray& vap, const SpatialPosition& pos, double t,
               double dt) { return p.value(va, vap, pos, t, dt); },
            py::arg("variable_array"),
            py::arg("variable_array_prev"),
            py::arg("pos"),
            py::arg("t"),
            py::arg("dt"),
            R"pbdoc(
             Evaluate the property value with previous time step data.

             Parameters:
                 variable_array: Current time step values.
                 variable_array_prev: Previous time step values.
                 pos: Spatial position.
                 t: Current time.
                 dt: Time step size.
             )pbdoc")

        .def(
            "dValue",
            [](const Property& p, const VariableArray& va,
               const VariableArray& vap, Variable var,
               const SpatialPosition& pos, double t, double dt)
            { return p.dValue(va, vap, var, pos, t, dt); },
            py::arg("variable_array"),
            py::arg("variable_array_prev"),
            py::arg("variable"),
            py::arg("pos"),
            py::arg("t"),
            py::arg("dt"),
            R"pbdoc(
             Evaluate the derivative of the property with respect to a variable,
             using previous time step data.

             Parameters:
                 variable_array: Current time step values.
                 variable_array_prev: Previous time step values.
                 variable: Variable to differentiate with respect to.
                 pos: Spatial position.
                 t: Current time.
                 dt: Time step size.
             )pbdoc")

        .def(
            "dValue",
            [](const Property& p, const VariableArray& va, Variable var,
               const SpatialPosition& pos, double t, double dt)
            { return p.dValue(va, var, pos, t, dt); },
            py::arg("variable_array"),
            py::arg("variable"),
            py::arg("pos"),
            py::arg("t"),
            py::arg("dt"),
            R"pbdoc(
             Evaluate the derivative of the property with respect to a variable.

             Parameters:
                 variable_array: Current time step values.
                 variable: Variable to differentiate with respect to.
                 pos: Spatial position.
                 t: Current time.
                 dt: Time step size.
             )pbdoc");
}

void bindConstant(py::module_& m)
{
    py::class_<Constant, Property>(m, "Constant", R"pbdoc(
        Constant property with a fixed value.
    )pbdoc")
        .def(py::init<std::string, PropertyDataType>(),
             py::arg("name"),
             py::arg("value"),
             R"pbdoc(
             Construct a Constant property.

             Parameters:
                 name: Name of the property.
                 value: Constant value.
             )pbdoc");
}

void bindLinear(py::module_& m)
{
    py::class_<Linear, Property>(m, "Linear", R"pbdoc(
        Linearly interpolated property dependent on independent variables.
    )pbdoc")
        .def(py::init(
                 [](std::string name, PropertyDataType reference_value,
                    std::vector<std::tuple<Variable, VariableType,
                                           VariableType>> const& ivs)
                 {
                     auto makeIV = [](std::tuple<Variable, VariableType,
                                                 VariableType> const& iv)
                     { return std::make_from_tuple<IndependentVariable>(iv); };

                     auto independent_variables =
                         ivs | ranges::views::transform(makeIV) |
                         ranges::to<std::vector>();
                     return std::make_unique<Linear>(name, reference_value,
                                                     independent_variables);
                 }),
             py::arg("name"),
             py::arg("reference_value"),
             py::arg("independent_variables"),
             R"pbdoc(
             Construct a Linear property.

             Parameters:
                 name: Name of the property.
                 reference_value: Base property value.
                 independent_variables: List of tuples (variable, reference condition, slope), each defining a linear dependency on a variable.
             )pbdoc");
}
}  // namespace

PYBIND11_MODULE(mpl, m)
{
#ifndef NDEBUG
    BaseLib::initOGSLogger("all");
#else   // NDEBUG
    BaseLib::initOGSLogger("info");
#endif  // NDEBUG

    m.attr("__name__") = "ogs.mpl";
    m.doc() = "pybind11 ogs MPL bindings";

    bindSpatialPosition(m);
    bindVariableArray(m);
    bindVariableEnum(m);
    bindProperty(m);
    bindConstant(m);
    bindLinear(m);
}
