/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MFrontGeneric.h"

namespace MaterialLib::Solids::MFront
{
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

    OGS_FATAL("Unknown variable type {}.", v);
}

int getEquivalentPlasticStrainOffset(mgis::behaviour::Behaviour const& b)
{
    return mgis::behaviour::contains(b.isvs, "EquivalentPlasticStrain")
               ? mgis::behaviour::getVariableOffset(
                     b.isvs, "EquivalentPlasticStrain", b.hypothesis)
               : -1;
}

template <int DisplacementDim>
OGSMFrontTangentOperatorData tangentOperatorDataMFrontToOGS(
    std::vector<double> const& mfront_data,
    std::optional<
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>> const& Q,
    mgis::behaviour::Behaviour const& behaviour)
{
    using VT = mgis::behaviour::Variable::Type;
    using KM = MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;
    using KV = MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;

    std::vector<double> ogs_data(mfront_data.size());

    std::size_t offset = 0;
    constexpr auto kv_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);

    for (auto const& [var1, var2] : behaviour.to_blocks)
    {
        std::size_t size;

        auto const* const d_in = mfront_data.data() + offset;
        auto* const d_out = ogs_data.data() + offset;

        if (var1.type == VT::SCALAR && var2.type == VT::SCALAR)
        {
            size = 1;
            *d_out = *d_in;
        }
        else if (var1.type == VT::SCALAR && var2.type == VT::STENSOR)
        {
            assert(getVariableSize(var2, behaviour.hypothesis) == kv_size);
            size = kv_size;

            if (Q)
            {
                OGS_FATAL(
                    "Coordinate frame rotation not yet implemented for "
                    "dScalar/dSTensor.");
            }

            Eigen::Map<KV>{d_out} = eigenSwap45View(Eigen::Map<const KV>{d_in});
        }
        else if (var1.type == VT::STENSOR && var2.type == VT::SCALAR)
        {
            assert(getVariableSize(var1, behaviour.hypothesis) == kv_size);
            size = kv_size;

            if (Q)
            {
                OGS_FATAL(
                    "Coordinate frame rotation not yet implemented for "
                    "dSTensor/dScalar.");
            }

            Eigen::Map<KV>{d_out} = eigenSwap45View(Eigen::Map<const KV>{d_in});
        }
        else if (var1.type == VT::STENSOR && var2.type == VT::STENSOR)
        {
            assert(getVariableSize(var1, behaviour.hypothesis) == kv_size);
            assert(getVariableSize(var2, behaviour.hypothesis) == kv_size);
            size = kv_size * kv_size;

            if (Q)
            {
                Eigen::Map<KM>{d_out} =
                    *Q * eigenSwap45View(Eigen::Map<const KM>{d_in}) *
                    Q->transpose();
            }
            else
            {
                Eigen::Map<KM>{d_out} =
                    eigenSwap45View(Eigen::Map<const KM>{d_in});
            }
        }
        else
        {
            OGS_FATAL("unsupported variable type combination");
        }

        offset += size;
    }

    return {std::move(ogs_data)};
}

template OGSMFrontTangentOperatorData tangentOperatorDataMFrontToOGS<2>(
    std::vector<double> const& mfront_data,
    std::optional<MathLib::KelvinVector::KelvinMatrixType<2>> const& Q,
    mgis::behaviour::Behaviour const& behaviour);
template OGSMFrontTangentOperatorData tangentOperatorDataMFrontToOGS<3>(
    std::vector<double> const& mfront_data,
    std::optional<MathLib::KelvinVector::KelvinMatrixType<3>> const& Q,
    mgis::behaviour::Behaviour const& behaviour);

}  // namespace MaterialLib::Solids::MFront
