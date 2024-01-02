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

#include <MGIS/Behaviour/Variable.hxx>

#include "MaterialLib/MPL/Utils/Tensor.h"
#include "MaterialLib/MPL/VariableType.h"
#include "MathLib/KelvinVector.h"

namespace MaterialLib::Solids::MFront
{
/**
 * Holds meta data describing a variable that can be passed from OGS to MFront.
 *
 * For the description of the parameter \c Derived please see Strain below.
 */
template <typename Derived>
struct Variable
{
    /// The number of components of the variable.
    template <int DisplacementDim>
    static constexpr std::size_t size()
    {
        return rows<DisplacementDim>() * cols<DisplacementDim>();
    }

    /// The number of rows of the variable.
    template <int DisplacementDim>
    static constexpr std::size_t rows()
    {
        using T = mgis::behaviour::Variable::Type;
        switch (Derived::type)
        {
            case T::SCALAR:
                return 1;
            case T::VECTOR:
                return DisplacementDim;
            case T::STENSOR:
                return MathLib::KelvinVector::kelvin_vector_dimensions(
                    DisplacementDim);
            case T::TENSOR:
                return MaterialPropertyLib::tensorSize(DisplacementDim);
        }
    }

    /// The number of columns of the variable.
    template <int DisplacementDim>
    static constexpr std::size_t cols()
    {
        using T = mgis::behaviour::Variable::Type;
        switch (Derived::type)
        {
            case T::SCALAR:
                return 1;
            case T::VECTOR:
                return 1;
            case T::STENSOR:
                return 1;
            case T::TENSOR:
                return 1;
        }
    }
};

/// Meta data for strain.
struct Strain : Variable<Strain>
{
    /// The name of the variable in MFront.
    constexpr static const char* name = "Strain";

    /// The type of the variable in MFront.
    constexpr static mgis::behaviour::Variable::Type type =
        mgis::behaviour::Variable::Type::STENSOR;

    /// The VariableArray entry that holds this variable in OGS.
    ///
    /// \note Currently we always pass strain via mechanical_strain.
    constexpr static auto mpl_var =
        &MaterialPropertyLib::VariableArray::mechanical_strain;
};

/// Instance that can be used for overload resolution/template type deduction.
static constexpr Strain strain;

/// Meta data for Green-Lagrange strain.
struct GreenLagrangeStrain : Variable<GreenLagrangeStrain>
{
    /// The name of the variable in MFront.
    constexpr static const char* name = "GreenLagrangeStrain";

    /// The type of the variable in MFront.
    constexpr static mgis::behaviour::Variable::Type type =
        mgis::behaviour::Variable::Type::STENSOR;

    /// The VariableArray entry that holds this variable in OGS.
    ///
    /// \note Currently we always pass strain via mechanical_strain.
    constexpr static auto mpl_var =
        &MaterialPropertyLib::VariableArray::mechanical_strain;
};

/// Instance that can be used for overload resolution/template type deduction.
static constexpr GreenLagrangeStrain green_lagrange_strain;

/// Meta data for deformation gradient.
struct DeformationGradient : Variable<DeformationGradient>
{
    /// The name of the variable in MFront.
    constexpr static const char* name = "DeformationGradient";

    /// The type of the variable in MFront.
    constexpr static mgis::behaviour::Variable::Type type =
        mgis::behaviour::Variable::Type::TENSOR;

    /// The VariableArray entry that holds this variable in OGS.
    constexpr static auto mpl_var =
        &MaterialPropertyLib::VariableArray::deformation_gradient;
};

/// Instance that can be used for overload resolution/template type deduction.
static constexpr DeformationGradient deformation_gradient;

struct LiquidPressure : Variable<LiquidPressure>
{
    constexpr static const char* name = "LiquidPressure";

    constexpr static mgis::behaviour::Variable::Type type =
        mgis::behaviour::Variable::Type::SCALAR;

    constexpr static auto mpl_var =
        &MaterialPropertyLib::VariableArray::liquid_phase_pressure;
};

static constexpr LiquidPressure liquid_pressure;

struct Stress : Variable<Stress>
{
    constexpr static const char* name = "Stress";

    constexpr static mgis::behaviour::Variable::Type type =
        mgis::behaviour::Variable::Type::STENSOR;

    constexpr static auto mpl_var = &MaterialPropertyLib::VariableArray::stress;
};

static constexpr Stress stress;

struct SecondPiolaKirchhoffStress : Variable<SecondPiolaKirchhoffStress>
{
    constexpr static const char* name = "SecondPiolaKirchhoffStress";

    constexpr static mgis::behaviour::Variable::Type type =
        mgis::behaviour::Variable::Type::STENSOR;

    constexpr static auto mpl_var = &MaterialPropertyLib::VariableArray::stress;
};

static constexpr SecondPiolaKirchhoffStress second_piola_kirchhoff_stress;

struct Saturation : Variable<Saturation>
{
    constexpr static const char* name = "Saturation";

    constexpr static mgis::behaviour::Variable::Type type =
        mgis::behaviour::Variable::Type::SCALAR;

    constexpr static auto mpl_var =
        &MaterialPropertyLib::VariableArray::liquid_saturation;
};

static constexpr Saturation saturation;

struct Temperature : Variable<Temperature>
{
    constexpr static const char* name = "Temperature";

    constexpr static mgis::behaviour::Variable::Type type =
        mgis::behaviour::Variable::Type::SCALAR;

    constexpr static auto mpl_var =
        &MaterialPropertyLib::VariableArray::temperature;
};

static constexpr Temperature temperature;
}  // namespace MaterialLib::Solids::MFront
