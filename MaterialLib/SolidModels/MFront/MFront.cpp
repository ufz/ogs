/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MFront.h"

#include <MGIS/Behaviour/Integrate.hxx>

#include "NumLib/Exceptions.h"

namespace
{
/// Converts between OGSes and MFront's Kelvin vector indices.
constexpr std::ptrdiff_t OGSToMFront(std::ptrdiff_t i)
{
    // MFront: 11 22 33 12 13 23
    // OGS:    11 22 33 12 23 13
    if (i < 4)
        return i;
    if (i == 4)
        return 5;
    return 4;
}
/// Converts between OGSes and MFront's Kelvin vector indices.
constexpr std::ptrdiff_t MFrontToOGS(std::ptrdiff_t i)
{
    // MFront: 11 22 33 12 13 23
    // OGS:    11 22 33 12 23 13
    return OGSToMFront(i);  // Same algorithm: indices 4 and 5 swapped.
}

/// Converts between OGSes and MFront's Kelvin vectors and matrices.
template <typename Derived>
typename Derived::PlainObject OGSToMFront(Eigen::DenseBase<Derived> const& m)
{
    static_assert(Derived::RowsAtCompileTime != Eigen::Dynamic, "Error");
    static_assert(Derived::ColsAtCompileTime != Eigen::Dynamic, "Error");

    typename Derived::PlainObject n;

    // optimal for row-major storage order
    for (std::ptrdiff_t r = 0; r < Eigen::DenseBase<Derived>::RowsAtCompileTime;
         ++r)
    {
        auto const R = OGSToMFront(r);
        for (std::ptrdiff_t c = 0;
             c < Eigen::DenseBase<Derived>::ColsAtCompileTime;
             ++c)
        {
            auto const C = OGSToMFront(c);
            n(R, C) = m(r, c);
        }
    }

    return n;
}

/// Converts between OGSes and MFront's Kelvin vectors and matrices.
template <typename Derived>
typename Derived::PlainObject MFrontToOGS(Eigen::DenseBase<Derived> const& m)
{
    static_assert(Derived::RowsAtCompileTime != Eigen::Dynamic, "Error");
    static_assert(Derived::ColsAtCompileTime != Eigen::Dynamic, "Error");

    typename Derived::PlainObject n;

    // optimal for row-major storage order
    for (std::ptrdiff_t r = 0; r < Eigen::DenseBase<Derived>::RowsAtCompileTime;
         ++r)
    {
        auto const R = MFrontToOGS(r);
        for (std::ptrdiff_t c = 0;
             c < Eigen::DenseBase<Derived>::ColsAtCompileTime;
             ++c)
        {
            auto const C = MFrontToOGS(c);
            n(R, C) = m(r, c);
        }
    }

    return n;
}

}  // namespace

namespace MaterialLib
{
namespace Solids
{
namespace MFront
{
const char* toString(mgis::behaviour::Behaviour::Kinematic kin)
{
    using K = mgis::behaviour::Behaviour::Kinematic;
    switch (kin)
    {
        case K::UNDEFINEDKINEMATIC:
            return "UNDEFINEDKINEMATIC";
        case K::SMALLSTRAINKINEMATIC:
            return "SMALLSTRAINKINEMATIC";
        case K::COHESIVEZONEKINEMATIC:
            return "COHESIVEZONEKINEMATIC";
        case K::FINITESTRAINKINEMATIC_F_CAUCHY:
            return "FINITESTRAINKINEMATIC_F_CAUCHY";
        case K::FINITESTRAINKINEMATIC_ETO_PK1:
            return "FINITESTRAINKINEMATIC_ETO_PK1";
    }

    OGS_FATAL("Unknown kinematic %d.", kin);
}

const char* toString(mgis::behaviour::Behaviour::Symmetry sym)
{
    using S = mgis::behaviour::Behaviour::Symmetry;
    switch (sym)
    {
        case S::ISOTROPIC:
            return "ISOTROPIC";
        case S::ORTHOTROPIC:
            return "ORTHOTROPIC";
    }

    OGS_FATAL("Unknown symmetry %d.", sym);
}
const char* btypeToString(int btype)
{
    using B = mgis::behaviour::Behaviour;
    if (btype == B::GENERALBEHAVIOUR)
        return "GENERALBEHAVIOUR";
    if (btype == B::STANDARDSTRAINBASEDBEHAVIOUR)
        return "STANDARDSTRAINBASEDBEHAVIOUR";
    if (btype == B::STANDARDFINITESTRAINBEHAVIOUR)
        return "STANDARDFINITESTRAINBEHAVIOUR";
    if (btype == B::COHESIVEZONEMODEL)
        return "COHESIVEZONEMODEL";

    OGS_FATAL("Unknown behaviour type %d.", btype);
}
const char* varTypeToString(int v)
{
    using V = mgis::behaviour::Variable;
    if (v == V::SCALAR)
        return "SCALAR";
    if (v == V::VECTOR)
        return "VECTOR";
    if (v == V::STENSOR)
        return "STENSOR";
    if (v == V::TENSOR)
        return "TENSOR";

    OGS_FATAL("Unknown variable type %d.", v);
}

template <int DisplacementDim>
MFront<DisplacementDim>::MFront(
    mgis::behaviour::Behaviour&& behaviour,
    std::vector<ParameterLib::Parameter<double> const*>&& material_properties)
    : _behaviour(std::move(behaviour)),
      _material_properties(std::move(material_properties))
{
    auto const hypothesis = behaviour.hypothesis;

    if (_behaviour.symmetry != mgis::behaviour::Behaviour::Symmetry::ISOTROPIC)
        OGS_FATAL(
            "The storage order of the stiffness matrix is not tested, yet. "
            "Thus, we cannot be sure if we compute the behaviour of "
            "anisotropic materials correctly. Therefore, currently only "
            "isotropic materials are allowed.");

    if (_behaviour.gradients.size() != 1)
        OGS_FATAL(
            "The behaviour must have exactly a single gradient as input.");

    if (_behaviour.gradients[0].name != "Strain")
        OGS_FATAL("The behaviour must be driven by strain.");

    if (_behaviour.gradients[0].type != mgis::behaviour::Variable::STENSOR)
        OGS_FATAL("Strain must be a symmetric tensor.");

    if (mgis::behaviour::getVariableSize(_behaviour.gradients[0], hypothesis) !=
        MFront<DisplacementDim>::KelvinVector::SizeAtCompileTime)
        OGS_FATAL("Strain must have %ld components instead of %lu.",
                  MFront<DisplacementDim>::KelvinVector::SizeAtCompileTime,
                  mgis::behaviour::getVariableSize(_behaviour.gradients[0],
                                                   hypothesis));

    if (_behaviour.thermodynamic_forces.size() != 1)
        OGS_FATAL(
            "The behaviour must compute exactly one thermodynamic force.");

    if (_behaviour.thermodynamic_forces[0].name != "Stress")
        OGS_FATAL("The behaviour must compute stress.");

    if (_behaviour.thermodynamic_forces[0].type !=
        mgis::behaviour::Variable::STENSOR)
        OGS_FATAL("Stress must be a symmetric tensor.");

    if (mgis::behaviour::getVariableSize(_behaviour.thermodynamic_forces[0],
                                         hypothesis) !=
        MFront<DisplacementDim>::KelvinVector::SizeAtCompileTime)
        OGS_FATAL("Stress must have %ld components instead of %lu.",
                  MFront<DisplacementDim>::KelvinVector::SizeAtCompileTime,
                  mgis::behaviour::getVariableSize(
                      _behaviour.thermodynamic_forces[0], hypothesis));

    if (!_behaviour.esvs.empty())
    {
        if (_behaviour.esvs[0].name != "Temperature")
        {
            OGS_FATAL(
                "Only temperature is supported as external state variable.");
        }

        if (mgis::behaviour::getVariableSize(_behaviour.esvs[0], hypothesis) !=
            1)
            OGS_FATAL(
                "Temperature must be a scalar instead of having %lu "
                "components.",
                mgis::behaviour::getVariableSize(
                    _behaviour.thermodynamic_forces[0], hypothesis));
    }

    if (_behaviour.mps.size() != _material_properties.size())
    {
        ERR("There are %d material properties in the loaded behaviour:",
            _behaviour.mps.size());
        for (auto const& mp : _behaviour.mps)
        {
            ERR("\t%s", mp.name.c_str());
        }
        OGS_FATAL("But the number of passed material properties is %d.",
                  _material_properties.size());
    }
}

template <int DisplacementDim>
std::unique_ptr<typename MechanicsBase<DisplacementDim>::MaterialStateVariables>
MFront<DisplacementDim>::createMaterialStateVariables() const
{
    return std::make_unique<MaterialStateVariables>(_behaviour);
}

template <int DisplacementDim>
std::optional<std::tuple<typename MFront<DisplacementDim>::KelvinVector,
                         std::unique_ptr<typename MechanicsBase<
                             DisplacementDim>::MaterialStateVariables>,
                         typename MFront<DisplacementDim>::KelvinMatrix>>
MFront<DisplacementDim>::integrateStress(
    double const t,
    ParameterLib::SpatialPosition const& x,
    double const dt,
    KelvinVector const& eps_prev,
    KelvinVector const& eps,
    KelvinVector const& sigma_prev,
    typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
        material_state_variables,
    double const T) const
{
    assert(
        dynamic_cast<MaterialStateVariables const*>(&material_state_variables));
    // New state, copy of current one, packed in unique_ptr for return.
    auto state = std::make_unique<MaterialStateVariables>(
        static_cast<MaterialStateVariables const&>(material_state_variables));
    auto& behaviour_data = state->_behaviour_data;

    // TODO add a test of material behaviour where the value of dt matters.
    behaviour_data.dt = dt;
    behaviour_data.rdt = 1.0;
    behaviour_data.K[0] = 4.0;  // if K[0] is greater than 3.5, the consistent
                                // tangent operator must be computed.

    // evaluate parameters at (t, x)
    {
        auto out = behaviour_data.s1.material_properties.begin();
        for (auto* param : _material_properties)
        {
            auto const& vals = (*param)(t, x);
            out = std::copy(vals.begin(), vals.end(), out);
        }
        assert(out == behaviour_data.s1.material_properties.end());
    }

    if (!behaviour_data.s1.external_state_variables.empty())
    {
        // assuming that there is only temperature
        behaviour_data.s1.external_state_variables[0] = T;
    }

    auto v = mgis::behaviour::make_view(behaviour_data);

    auto const eps_prev_MFront = OGSToMFront(eps_prev);
    for (auto i = 0; i < KelvinVector::SizeAtCompileTime; ++i)
    {
        v.s0.gradients[i] = eps_prev_MFront[i];
    }

    auto const eps_MFront = OGSToMFront(eps);
    for (auto i = 0; i < KelvinVector::SizeAtCompileTime; ++i)
    {
        v.s1.gradients[i] = eps_MFront[i];
    }

    auto const sigma_prev_MFront = OGSToMFront(sigma_prev);
    for (auto i = 0; i < KelvinVector::SizeAtCompileTime; ++i)
    {
        v.s0.thermodynamic_forces[i] = sigma_prev_MFront[i];
        v.s1.thermodynamic_forces[i] = sigma_prev_MFront[i];
    }

    auto const status = mgis::behaviour::integrate(v, _behaviour);
    if (status != 1)
    {
        throw NumLib::AssemblyException("MFront: integration failed with status"
                + std::to_string(status) + ".");
    }

    KelvinVector sigma;
    for (auto i = 0; i < KelvinVector::SizeAtCompileTime; ++i)
    {
        sigma[i] = behaviour_data.s1.thermodynamic_forces[i];
    }
    sigma = MFrontToOGS(sigma);

    // TODO row- vs. column-major storage order. This should only matter for
    // anisotropic materials.
    if (behaviour_data.K.size() !=
        KelvinMatrix::RowsAtCompileTime * KelvinMatrix::ColsAtCompileTime)
        OGS_FATAL("Stiffness matrix has wrong size.");

    KelvinMatrix C =
        MFrontToOGS(Eigen::Map<KelvinMatrix>(behaviour_data.K.data()));

    return std::make_optional(
        std::make_tuple<typename MFront<DisplacementDim>::KelvinVector,
                        std::unique_ptr<typename MechanicsBase<
                            DisplacementDim>::MaterialStateVariables>,
                        typename MFront<DisplacementDim>::KelvinMatrix>(
            std::move(sigma), std::move(state), std::move(C)));
}

template <int DisplacementDim>
std::vector<typename MechanicsBase<DisplacementDim>::InternalVariable>
MFront<DisplacementDim>::getInternalVariables() const
{
    std::vector<typename MechanicsBase<DisplacementDim>::InternalVariable>
        internal_variables;

    for (auto const& iv : _behaviour.isvs)
    {
        auto const name = iv.name;
        auto const offset = mgis::behaviour::getVariableOffset(
            _behaviour.isvs, name, _behaviour.hypothesis);
        auto const size =
            mgis::behaviour::getVariableSize(iv, _behaviour.hypothesis);

        typename MechanicsBase<DisplacementDim>::InternalVariable new_variable{
            name, static_cast<unsigned>(size),
            [offset, size](
                typename MechanicsBase<
                    DisplacementDim>::MaterialStateVariables const& state,
                std::vector<double>& cache) -> std::vector<double> const& {
                assert(dynamic_cast<MaterialStateVariables const*>(&state) !=
                       nullptr);
                auto const& internal_state_variables =
                    static_cast<MaterialStateVariables const&>(state)
                        ._behaviour_data.s1.internal_state_variables;

                cache.resize(size);
                std::copy_n(internal_state_variables.data() + offset,
                            size,
                            begin(cache));
                return cache;
            }};
        internal_variables.push_back(new_variable);
    }

    return internal_variables;
}

template <int DisplacementDim>
double MFront<DisplacementDim>::computeFreeEnergyDensity(
    double const /*t*/,
    ParameterLib::SpatialPosition const& /*x*/,
    double const /*dt*/,
    KelvinVector const& /*eps*/,
    KelvinVector const& /*sigma*/,
    typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
    /*material_state_variables*/) const
{
    // TODO implement
    return std::numeric_limits<double>::quiet_NaN();
}

template class MFront<2>;
template class MFront<3>;

}  // namespace MFront
}  // namespace Solids
}  // namespace MaterialLib
