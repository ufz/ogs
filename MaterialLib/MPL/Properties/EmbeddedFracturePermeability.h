/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include "MaterialLib/MPL/Property.h"
#include "MathLib/KelvinVector.h"

namespace MaterialPropertyLib
{
/**
 * \class EmbeddedFracturePermeability
 * \brief Permeability model as proposed by Olivella&Alonso
 * \details This property must be a medium property, it
 * computes the permeability in dependence of the strain
 */
template <int DisplacementDim>
class EmbeddedFracturePermeability final : public Property
{
private:
    Medium* _medium = nullptr;
    Eigen::Matrix<double, 3, 1> const _n;
    bool const _n_const;
    double const _k;
    double const _b0;
    double const _a;
    double const _e0;

public:
    EmbeddedFracturePermeability(
        std::string name,
        Eigen::Matrix<double, 3, 1> const fracture_normal,
        bool const fracture_normal_is_constant,
        double const intrinsic_permeability,
        double const initial_aperture,
        double const mean_fracture_distance,
        double const threshold_strain);

    static int const KelvinVectorSize =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    using SymmetricTensor = Eigen::Matrix<double, KelvinVectorSize, 1>;

    void checkScale() const override;

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const primary_variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;
};

}  // namespace MaterialPropertyLib
