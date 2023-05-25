/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <MGIS/Behaviour/Behaviour.hxx>
#include <MGIS/Behaviour/BehaviourData.hxx>
#include <MGIS/Behaviour/Integrate.hxx>
#include <boost/mp11.hpp>

#include "MaterialLib/MPL/Utils/GetSymmetricTensor.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "NumLib/Exceptions.h"
#include "ParameterLib/Parameter.h"
#include "TangentOperatorBlocksView.h"
#include "ThermodynamicForcesView.h"

namespace MaterialLib::Solids::MFront
{
namespace MPL = MaterialPropertyLib;

/// Converts between OGSes and MFront's Kelvin vector indices.
std::ptrdiff_t OGSToMFront(std::ptrdiff_t i);

/// Converts between OGSes and MFront's Kelvin vector indices.
std::ptrdiff_t MFrontToOGS(std::ptrdiff_t i);

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

const char* varTypeToString(int v);

int getEquivalentPlasticStrainOffset(mgis::behaviour::Behaviour const& b);

/** Transforms MFront's to OGS's tangent operator data.
 *
 * Essentially swaps some off-diagonal components of symmetric tensors and
 * applies the given rotation of the frame of reference \c Q.
 */
template <int DisplacementDim>
OGSMFrontTangentOperatorData tangentOperatorDataMFrontToOGS(
    std::vector<double> const& mfront_data,
    std::optional<
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>> const& Q,
    mgis::behaviour::Behaviour const& behaviour);

extern template OGSMFrontTangentOperatorData tangentOperatorDataMFrontToOGS<2>(
    std::vector<double> const& mfront_data,
    std::optional<MathLib::KelvinVector::KelvinMatrixType<2>> const& Q,
    mgis::behaviour::Behaviour const& behaviour);
extern template OGSMFrontTangentOperatorData tangentOperatorDataMFrontToOGS<3>(
    std::vector<double> const& mfront_data,
    std::optional<MathLib::KelvinVector::KelvinMatrixType<3>> const& Q,
    mgis::behaviour::Behaviour const& behaviour);

// TODO template parameter only because of base class.
template <int DisplacementDim>
struct MaterialStateVariablesMFront
    : public MechanicsBase<DisplacementDim>::MaterialStateVariables
{
    explicit MaterialStateVariablesMFront(
        int const equivalent_plastic_strain_offset,
        mgis::behaviour::Behaviour const& b)
        : equivalent_plastic_strain_offset_(equivalent_plastic_strain_offset),
          _behaviour_data{b}
    {
    }

    MaterialStateVariablesMFront(
        MaterialStateVariablesMFront<DisplacementDim> const&) = default;
    MaterialStateVariablesMFront(
        MaterialStateVariablesMFront<DisplacementDim>&&) = delete;

    void pushBackState() override { mgis::behaviour::update(_behaviour_data); }

    int const equivalent_plastic_strain_offset_;
    mgis::behaviour::BehaviourData _behaviour_data;

    double getEquivalentPlasticStrain() const override
    {
        if (equivalent_plastic_strain_offset_ >= 0)
        {
            return _behaviour_data.s1
                .internal_state_variables[static_cast<mgis::size_type>(
                    equivalent_plastic_strain_offset_)];
        }

        return 0.0;
    }
};

namespace detail
{
/// A struct mapping a MGIS variable type to its corresponding OGS type.
template <int DisplacementDim, mgis::behaviour::Variable::Type MFrontType>
struct MapToMPLType;

template <int DisplacementDim>
struct MapToMPLType<DisplacementDim, mgis::behaviour::Variable::Type::STENSOR>
{
    using type = MaterialPropertyLib::SymmetricTensor<DisplacementDim>;
};

template <int DisplacementDim>
struct MapToMPLType<DisplacementDim, mgis::behaviour::Variable::Type::SCALAR>
{
    using type = double;
};

template <int DisplacementDim, mgis::behaviour::Variable::Type MFrontType>
using MapToMPLType_t = typename MapToMPLType<DisplacementDim, MFrontType>::type;

/// A helper struct filling MFront's gradient data (or thermodynamic forces)
/// with values taken from a MaterialPropertyLib::VariableArray.
template <int DisplacementDim>
struct SetGradient
{
    MaterialPropertyLib::VariableArray const& variable_array;
    std::optional<
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>> const& Q;
    double* target;

    template <typename Grad>
    void operator()(Grad)
    {
        auto constexpr num_comp = Grad::template size<DisplacementDim>();

        if constexpr (Grad::type == mgis::behaviour::Variable::Type::SCALAR)
        {
            *target = variable_array.*Grad::mpl_var;
        }
        else if constexpr (Grad::type ==
                           mgis::behaviour::Variable::Type::STENSOR)
        {
            using MPLType = MapToMPLType_t<DisplacementDim, Grad::type>;
            auto const& grad_ogs =
                std::get<MPLType>(variable_array.*Grad::mpl_var);

            auto const grad_mfront = Q ? OGSToMFront(Q->transpose() * grad_ogs)
                                       : OGSToMFront(grad_ogs);
            std::copy_n(grad_mfront.data(), num_comp, target);
        }
        else
        {
            OGS_FATAL("Unsupported gradient type {}.",
                      varTypeToString(Grad::type));
        }

        target += num_comp;
    }
};
}  // namespace detail

/// Uses a material model provided by MFront (via MFront's generic interface and
/// the MGIS library).
template <int DisplacementDim,
          typename Gradients,
          typename TDynForces,
          typename ExtStateVars>
class MFrontGeneric
{
    static_assert(boost::mp11::mp_is_set<Gradients>::value);
    static_assert(boost::mp11::mp_is_set<TDynForces>::value);
    static_assert(boost::mp11::mp_is_set<ExtStateVars>::value);

    using GradientsAndExtStateVars =
        boost::mp11::mp_append<Gradients, ExtStateVars>;

    static_assert(boost::mp11::mp_is_set<GradientsAndExtStateVars>::value);

public:
    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

    MFrontGeneric(mgis::behaviour::Behaviour&& behaviour,
                  std::vector<ParameterLib::Parameter<double> const*>&&
                      material_properties,
                  std::optional<ParameterLib::CoordinateSystem> const&
                      local_coordinate_system)
        : _behaviour(std::move(behaviour)),
          equivalent_plastic_strain_offset_(
              getEquivalentPlasticStrainOffset(_behaviour)),
          _material_properties(std::move(material_properties)),
          _local_coordinate_system(local_coordinate_system
                                       ? &local_coordinate_system.value()
                                       : nullptr)
    {
        {
            auto check_gradient = [&grads = _behaviour.gradients,
                                   hyp = _behaviour.hypothesis,
                                   i = 0]<typename Grad>(Grad) mutable
            {
                // TODO allow reordering of gradients and thermodynamic forces?
                if (grads[i].name != Grad::name)
                {
                    OGS_FATAL("The behaviour's {}th driver must be {}.", i,
                              Grad::name);
                }

                if (grads[i].type != Grad::type)
                {
                    OGS_FATAL(
                        "The behaviour's {}th driver ({}) must be of type {}.",
                        i, grads[i].name, varTypeToString(Grad::type));
                }

                if (mgis::behaviour::getVariableSize(grads[i], hyp) !=
                    Grad::template size<DisplacementDim>())
                {
                    OGS_FATAL(
                        "The behaviour's {}th driver ({}) must have size {} "
                        "instead of {}.",
                        i, grads[i].name,
                        Grad::template size<DisplacementDim>(),
                        mgis::behaviour::getVariableSize(grads[i], hyp));
                }

                i++;
            };

            if (_behaviour.gradients.size() !=
                boost::mp11::mp_size<Gradients>::value)
                OGS_FATAL(
                    "The behaviour must have exactly {} gradients as input.",
                    boost::mp11::mp_size<Gradients>::value);

            boost::mp11::mp_for_each<Gradients>(check_gradient);
        }

        {
            auto check_tdyn_force = [&tdfs = _behaviour.thermodynamic_forces,
                                     hyp = _behaviour.hypothesis,
                                     i = 0]<typename TDF>(TDF) mutable
            {
                if (tdfs[i].name != TDF::name)
                {
                    OGS_FATAL(
                        "The behaviour's {}th thermodynamic force must be {}.",
                        i, TDF::name);
                }

                if (tdfs[i].type != TDF::type)
                {
                    OGS_FATAL(
                        "The behaviour's {}th thermodynamic force ({}) must be "
                        "of "
                        "type {}.",
                        i, tdfs[i].name, varTypeToString(TDF::type));
                }

                if (mgis::behaviour::getVariableSize(tdfs[i], hyp) !=
                    TDF::template size<DisplacementDim>())
                {
                    OGS_FATAL(
                        "The behaviour's {}th thermodynamic force ({}) must "
                        "have "
                        "size {} instead of {}.",
                        i, tdfs[i].name, TDF::template size<DisplacementDim>(),
                        mgis::behaviour::getVariableSize(tdfs[i], hyp));
                }

                i++;
            };

            if (_behaviour.thermodynamic_forces.size() !=
                boost::mp11::mp_size<TDynForces>::value)
                OGS_FATAL(
                    "The behaviour must compute exactly {} thermodynamic "
                    "forces.",
                    boost::mp11::mp_size<TDynForces>::value);

            boost::mp11::mp_for_each<TDynForces>(check_tdyn_force);
        }

        auto const hypothesis = _behaviour.hypothesis;

        if (!_behaviour.esvs.empty())
        {
            if (_behaviour.esvs[0].name != "Temperature")
            {
                OGS_FATAL(
                    "Only temperature is supported as external state "
                    "variable.");
            }

            if (mgis::behaviour::getVariableSize(_behaviour.esvs[0],
                                                 hypothesis) != 1)
                OGS_FATAL(
                    "Temperature must be a scalar instead of having {:d} "
                    "components.",
                    mgis::behaviour::getVariableSize(
                        _behaviour.thermodynamic_forces[0], hypothesis));
        }

        if (_behaviour.mps.size() != _material_properties.size())
        {
            ERR("There are {:d} material properties in the loaded behaviour:",
                _behaviour.mps.size());
            for (auto const& mp : _behaviour.mps)
            {
                ERR("\t{:s}", mp.name);
            }
            OGS_FATAL("But the number of passed material properties is {:d}.",
                      _material_properties.size());
        }
    }

    std::unique_ptr<
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables>
    createMaterialStateVariables() const
    {
        return std::make_unique<MaterialStateVariablesMFront<DisplacementDim>>(
            equivalent_plastic_strain_offset_, _behaviour);
    }

    std::optional<std::tuple<OGSMFrontThermodynamicForcesData,
                             std::unique_ptr<typename MechanicsBase<
                                 DisplacementDim>::MaterialStateVariables>,
                             OGSMFrontTangentOperatorData>>
    integrateStress(
        MaterialPropertyLib::VariableArray const& variable_array_prev,
        MaterialPropertyLib::VariableArray const& variable_array,
        double const t,
        ParameterLib::SpatialPosition const& x,
        double const dt,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables) const
    {
        using namespace MathLib::KelvinVector;

        assert(
            dynamic_cast<MaterialStateVariablesMFront<DisplacementDim> const*>(
                &material_state_variables));
        // New state, copy of current one, packed in unique_ptr for return.
        auto state = std::make_unique<
            MaterialStateVariablesMFront<DisplacementDim>>(
            static_cast<MaterialStateVariablesMFront<DisplacementDim> const&>(
                material_state_variables));
        auto& behaviour_data = state->_behaviour_data;

        behaviour_data.dt = dt;
        behaviour_data.rdt = 1.0;
        behaviour_data.K[0] =
            4.0;  // if K[0] is greater than 3.5, the consistent
                  // tangent operator must be computed.

        // evaluate parameters at (t, x)
        {
            {
                auto out = behaviour_data.s0.material_properties.begin();
                for (auto* param : _material_properties)
                {
                    auto const& vals = (*param)(t - dt, x);
                    out = std::copy(vals.begin(), vals.end(), out);
                }
                assert(out == behaviour_data.s0.material_properties.end());
            }

            {
                auto out = behaviour_data.s1.material_properties.begin();
                for (auto* param : _material_properties)
                {
                    auto const& vals = (*param)(t, x);
                    out = std::copy(vals.begin(), vals.end(), out);
                }
                assert(out == behaviour_data.s1.material_properties.end());
            }
        }

        // TODO unify with gradient handling? Make gradient and external state
        // var input both optional?
        if (!behaviour_data.s0.external_state_variables.empty())
        {
            // assuming that there is only temperature
            // NOTE the temperature can be NaN, e.g., if OGS's process does not
            // have a temperature defined
            behaviour_data.s0.external_state_variables[0] =
                variable_array_prev.temperature;
        }

        if (!behaviour_data.s1.external_state_variables.empty())
        {
            // assuming that there is only temperature
            // NOTE the temperature can be NaN, e.g., if OGS's process does not
            // have a temperature defined
            behaviour_data.s1.external_state_variables[0] =
                variable_array.temperature;
        }

        // rotation tensor
        std::optional<KelvinMatrixType<DisplacementDim>> Q;
        if (_local_coordinate_system)
        {
            Q = fourthOrderRotationMatrix(
                _local_coordinate_system->transformation<DisplacementDim>(x));
        }

        boost::mp11::mp_for_each<Gradients>(
            detail::SetGradient<DisplacementDim>{
                variable_array_prev, Q, behaviour_data.s0.gradients.data()});
        boost::mp11::mp_for_each<Gradients>(
            detail::SetGradient<DisplacementDim>{
                variable_array, Q, behaviour_data.s1.gradients.data()});

        // previous and current state of thermodynamic forces are both set from
        // variable_array_prev. TODO optimization potential compute Q * grad_ogs
        // only once for both cases.
        boost::mp11::mp_for_each<TDynForces>(
            detail::SetGradient<DisplacementDim>{
                variable_array_prev, Q,
                behaviour_data.s0.thermodynamic_forces.data()});
        boost::mp11::mp_for_each<TDynForces>(
            detail::SetGradient<DisplacementDim>{
                variable_array_prev, Q,
                behaviour_data.s1.thermodynamic_forces.data()});

        auto v = mgis::behaviour::make_view(behaviour_data);
        auto const status = mgis::behaviour::integrate(v, _behaviour);
        if (status != 1)
        {
            throw NumLib::AssemblyException(
                "MFront: integration failed with status " +
                std::to_string(status) + ".");
        }

        OGSMFrontThermodynamicForcesData tdyn_forces_data;
        tdyn_forces_data.data.resize(  // TODO data stored on heap but size is
                                       // compile-time constant
            behaviour_data.s1.thermodynamic_forces.size());

        boost::mp11::mp_for_each<TDynForces>(
            [&out_data = tdyn_forces_data,
             &in_data = behaviour_data.s1.thermodynamic_forces,
             &Q]<typename TDF>(TDF tdf)
            {
                OGSMFrontThermodynamicForcesView<DisplacementDim, TDynForces>
                    view;

                if constexpr (TDF::type ==
                              mgis::behaviour::Variable::Type::STENSOR)
                {
                    if (Q)
                    {
                        view.block(tdf, out_data) =
                            *Q * MFrontToOGS(view.block(tdf, in_data));
                    }
                    else
                    {
                        view.block(tdf, out_data) =
                            MFrontToOGS(view.block(tdf, in_data));
                    }
                }
                else if constexpr (TDF::type ==
                                   mgis::behaviour::Variable::Type::SCALAR)
                {
                    view.block(tdf, out_data) = view.block(tdf, in_data);
                }
                else
                {
                    OGS_FATAL("Not yet implemented.");
                }
            });

        return std::make_optional(
            std::make_tuple<OGSMFrontThermodynamicForcesData,
                            std::unique_ptr<typename MechanicsBase<
                                DisplacementDim>::MaterialStateVariables>,
                            OGSMFrontTangentOperatorData>(
                std::move(tdyn_forces_data),
                std::move(state),
                tangentOperatorDataMFrontToOGS<DisplacementDim>(
                    behaviour_data.K, Q, _behaviour)));
    }

    using InternalVariable =
        typename MechanicsBase<DisplacementDim>::InternalVariable;

    std::vector<InternalVariable> getInternalVariables() const
    {
        std::vector<InternalVariable> internal_variables;

        for (auto const& iv : _behaviour.isvs)
        {
            auto const name = iv.name;
            auto const offset = mgis::behaviour::getVariableOffset(
                _behaviour.isvs, name, _behaviour.hypothesis);
            auto const size =
                mgis::behaviour::getVariableSize(iv, _behaviour.hypothesis);

            // TODO (naumov): For orthotropic materials the internal variables
            // should be rotated to the global coordinate system before output.
            // MFront stores the variables in local coordinate system.
            // The `size` variable could be used to find out the type of
            // variable.
            InternalVariable new_variable{
                name, static_cast<int>(size),
                [offset, size](
                    typename MechanicsBase<
                        DisplacementDim>::MaterialStateVariables const& state,
                    std::vector<double>& cache) -> std::vector<double> const&
                {
                    assert(dynamic_cast<MaterialStateVariablesMFront<
                               DisplacementDim> const*>(&state) != nullptr);
                    auto const& internal_state_variables =
                        static_cast<MaterialStateVariablesMFront<
                            DisplacementDim> const&>(state)
                            ._behaviour_data.s1.internal_state_variables;

                    cache.resize(size);
                    std::copy_n(internal_state_variables.data() + offset,
                                size,
                                begin(cache));
                    return cache;
                },
                [offset, size](typename MechanicsBase<
                               DisplacementDim>::MaterialStateVariables& state)
                    -> std::span<double>
                {
                    assert(dynamic_cast<MaterialStateVariablesMFront<
                               DisplacementDim> const*>(&state) != nullptr);
                    auto& internal_state_variables =
                        static_cast<
                            MaterialStateVariablesMFront<DisplacementDim>&>(
                            state)
                            ._behaviour_data.s1.internal_state_variables;

                    return {internal_state_variables.data() + offset, size};
                }};
            internal_variables.push_back(new_variable);
        }

        return internal_variables;
    }

    OGSMFrontTangentOperatorBlocksView<DisplacementDim,
                                       Gradients,
                                       TDynForces,
                                       ExtStateVars>
    createTangentOperatorBlocksView() const
    {
        return OGSMFrontTangentOperatorBlocksView<DisplacementDim,
                                                  Gradients,
                                                  TDynForces,
                                                  ExtStateVars>{
            _behaviour.to_blocks};
    }

    OGSMFrontThermodynamicForcesView<DisplacementDim, TDynForces>
    createThermodynamicForcesView() const
    {
        return OGSMFrontThermodynamicForcesView<DisplacementDim, TDynForces>{};
    }

    double getBulkModulus(double const /*t*/,
                          ParameterLib::SpatialPosition const& /*x*/,
                          KelvinMatrix const* const C) const
    {
        if (C == nullptr)
        {
            OGS_FATAL(
                "MFront::getBulkModulus() requires the tangent stiffness C "
                "input "
                "argument to be valid.");
        }
        auto const& identity2 = MathLib::KelvinVector::Invariants<
            MathLib::KelvinVector::kelvin_vector_dimensions(
                DisplacementDim)>::identity2;
        return 1. / 9. * identity2.transpose() * *C * identity2;
    }

    double computeFreeEnergyDensity(
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

private:
    mgis::behaviour::Behaviour _behaviour;
    int const equivalent_plastic_strain_offset_;
    std::vector<ParameterLib::Parameter<double> const*> _material_properties;
    ParameterLib::CoordinateSystem const* const _local_coordinate_system;
};
}  // namespace MaterialLib::Solids::MFront
