/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on November 10, 2020, 8:49 AM
 */

#pragma once

#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/VariableType.h"

namespace ParameterLib
{
struct CoordinateSystem;
template <typename T>
struct Parameter;
}  // namespace ParameterLib

namespace MaterialPropertyLib
{
/**
 * \brief A strain dependent intrinsic permeability model.
 *
 *  The model was proposed
 *  in \cite xu2011simulation and it was further investigated
 *  in \cite xu2013coupled .
 *
 *   The model takes the form of
 *  \f[ \mathbf{k} =f(\epsilon_v) e^{b_1  {\bar\epsilon}^p}\mathbf{k}_0
 *  \f]
 *  with
 *  \f[ f(\epsilon_v)=
 *      \begin{cases}
 *       10^{b_2  \epsilon_v}, & \epsilon_v <=0\\
 *      10^{b_3  \epsilon_v}, & \epsilon_v >0
 *      \end{cases}
 *  \f]
 *   where
 *   <table>
 *   <tr><td>\f$ \epsilon_v \f$  <td> the volumetric strain,
 *   <tr><td> \f$ {\bar\epsilon}^p\f$ <td> the equivalent plastic
 *    strain,
 *   <tr><td>\f$\mathbf{k}_0\f$  <td> the initial intrinsic permeability,
 *   <tr><td>\f$b_1,\,b_2,\,b_3\f$  <td> the three parameters.
 * </table>
 *
 *  * Note: In \cite xu2011simulation  and \cite xu2013coupled ,
 * from the point of view of experiment of permeability change, the symbols of
 * \f$ \Delta \epsilon_v \f$ and \f$\Delta {\bar\epsilon}^p\f$  are used for \f$
 * \epsilon_v \f$ and \f$
 * {\bar\epsilon}^p\f$, respectively. That means that the symbols of
 * \f$ \Delta \epsilon_v \f$ and \f$\Delta {\bar\epsilon}^p\f$
 *  refer to the increment of the volumetric strain and the increment of the
 * equivalent plastic strain, respectively, from the beginning of the
 * experiment.
 */
template <int DisplacementDim>
class StrainDependentPermeability final : public Property
{
public:
    StrainDependentPermeability(
        std::string name, ParameterLib::Parameter<double> const& k0,
        double const b1, double const b2, double const b3,
        double const minimum_permeability, double const maximum_permeability,
        ParameterLib::CoordinateSystem const* const local_coordinate_system);

    void checkScale() const override;

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;

private:
    /// Initial intrinsic permeability.
    ParameterLib::Parameter<double> const& k0_;
    double const b1_;
    double const b2_;
    double const b3_;
    double const minimum_permeability_;
    double const maximum_permeability_;
    ParameterLib::CoordinateSystem const* const local_coordinate_system_;
};

extern template class StrainDependentPermeability<2>;
extern template class StrainDependentPermeability<3>;

}  // namespace MaterialPropertyLib
