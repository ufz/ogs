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
#include "ParameterLib/Parameter.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;
/**
 * \class OrthotropicEmbeddedFracturePermeability
 * \brief Extended Permeability model based on Olivella&Alonso
 * \details This property must be a medium property, it
 * computes the permeability in dependence of the strain
 *
 * The base model was proposed
 * in \cite alonso2006mechanisms and it was further investigated
 * in \cite olivella2008gas . This extended Version features three
 * orthotropic fracture planes.
 *
 * The model takes the form of
 * \f[ \mathbf{k} = k_\text{m} \mathbf{I} + \sum \limits_{i=1}^3
 * \frac{b_i}{a_i} \left( \frac{b_i^2}{12} - k_\text{m} \right) \left(
 * \mathbf{I} - \mathbf{M}_i \right) \f]
 * with
 * \f[ \mathbf{M}_i = \vec{n}_i \otimes \vec{n}_i \f]
 * and
 * \f[ b_i = b_{i0} + \Delta b_i \\
 * \Delta b_i = a_i \langle \mathbf{\epsilon} : \mathbf{M}_i -
 * \varepsilon_{0i} \rangle
 * \f]
 * where
 * <table>
 * <tr><td>\f$ k_\text{m} \f$  <td> permeability of undisturbed material
 * <tr><td>\f$ b_i \f$  <td> fracture aperture of each fracture plane
 * <tr><td>\f$ b_{i0} \f$  <td> initial aperture of each fracture plane
 * <tr><td>\f$ a_i \f$  <td> mean fracture distance of each fracture plane
 * <tr><td>\f$ \vec{n}_i \f$  <td> fracture normal vector of each fracture
 * plane
 * <tr><td>\f$ \varepsilon_{i0} \f$  <td> threshold strain of each fracture
 * plane
 * </table>
 */

template <int DisplacementDim>
class OrthotropicEmbeddedFracturePermeability final : public Property
{
private:
    Medium* _medium = nullptr;
    std::vector<double> const _a;
    std::vector<double> const _e0;
    Eigen::Matrix<double, 3, 3> const _n;
    ParameterLib::Parameter<double> const& _k;
    ParameterLib::Parameter<double> const& _phi_xy;
    ParameterLib::Parameter<double> const& _phi_yz;

public:
    OrthotropicEmbeddedFracturePermeability(
        std::string name,
        std::vector<double> const& mean_fracture_distances,
        std::vector<double> const& threshold_strains,
        Eigen::Matrix<double, 3, 3> const fracture_normals,
        ParameterLib::Parameter<double> const& intrinsic_permeability,
        ParameterLib::Parameter<double> const& fracture_rotation_xy,
        ParameterLib::Parameter<double> const& fracture_rotation_yz);

    using SymmetricTensor = Eigen::Matrix<
        double,
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim), 1>;

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
