/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   LiquidFlowMaterialProperties.h
 *
 * Created on August 18, 2016, 11:03 AM
 */

#ifndef OGS_LIQUIDFLOWMATERIALPROPERTIES_H
#define OGS_LIQUIDFLOWMATERIALPROPERTIES_H

#include <memory>

#include "MaterialLib/Fluid/FluidPropertyHeaders.h"
#include "MaterialLib/PorousMedium/PorousPropertyHeaders.h"

namespace MaterialLib
{
namespace PorousMedium
{
class Porosity;
class Storage;
}
}

namespace MeshLib
{
template <typename PROP_VAL_TYPE> class PropertyVector;
}

namespace ProcessLib
{
template <typename T> struct Parameter;
class SpatialPosition;

namespace LiquidFlow
{
class LiquidFlowMaterialProperties
{
public:
    typedef MaterialLib::Fluid::FluidProperty::ArrayType ArrayType;

    LiquidFlowMaterialProperties(
            BaseLib::ConfigTree const& config,
            MeshLib::PropertyVector<int> const& material_ids,
            Parameter<double> const& intrinsic_permeability_data,
            Parameter<double> const& porosity_data,
            Parameter<double> const& storage_data);

    /**
     * \brief Compute the coefficient of the mass term by
     *      \f[
     *           n \frac{partial \rho_l}{\partial p} + \beta_s
     *      \f]
     *     where \f$n\f$ is the porosity, \f$rho_l\f$ is the liquid density,
     *     \f$bata_s\f$ is the storage.
     * \param t                  Time.
     * \param pos                Position of element.
     * \param dim                Dimension of space.
     * \param porosity_variable  The first variable for porosity model, and it
     *                           passes a double type value that could be
     *                           saturation, and invariant of stress or strain.
     * \param storage_variable   Variable for storage model.
     * \param p                  Pressure value
     * \param T                  Temperature value
     */
    double getMassCoefficient(const double t, const SpatialPosition& pos,
                              const double p, const double T,
                              const double porosity_variable,
                              const double storage_variable) const;

    Eigen::MatrixXd const& getPermeability(const double t,
                                            const SpatialPosition& pos,
                                            const int dim) const;

    double getLiquidDensity(const double p, const double T) const;

    double getViscosity(const double p, const double T) const;

private:
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> _liquid_density;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> _viscosity;

    /** Use porous medium models for different material zones.
     *  Material IDs must be given as mesh element properties.
     */
    MeshLib::PropertyVector<int> const& _material_ids;
    std::vector<Eigen::MatrixXd> _intrinsic_permeability_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>> _porosity_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Storage>> _storage_models;

    /// Use data for porous medium properties.
    Parameter<double> const& _intrinsic_permeability_data;
    Parameter<double> const& _porosity_data;
    Parameter<double> const& _storage_data;
};

}  // end of namespace
}  // end of namespace
#endif /* LIQUIDFLOWMATERIALPROPERTIES_H */
