/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 10, 2019
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BaseLib/ConfigTree.h"
#include "Constant.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Constant> createConstant(BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type", "Constant");
    DBUG("Create Constant property");
    std::vector<double> const values =
        //! \ogs_file_param{properties__property__Constant__value}
        config.getConfigParameter<std::vector<double>>("value");

    switch (values.size())
    {
        case 1:
        {
            // scalar
            PropertyDataType property_value = values[0];
            return std::make_unique<Constant>(property_value);
        }
        case 2:
        {
            // Pair
            PropertyDataType property_value = Pair{values[0], values[1]};
            return std::make_unique<Constant>(property_value);
        }
        case 3:
        {
            // Vector
            PropertyDataType property_value =
                Vector{values[0], values[1], values[2]};
            return std::make_unique<Constant>(property_value);
        }
        case 4:
        {
            // Tensor
            PropertyDataType property_value =
                Tensor2d{values[0], values[1], values[2], values[3]};
            return std::make_unique<Constant>(property_value);
        }
        case 6:
        {
            // Symmetric Tensor - xx, yy, zz, xy, xz, yz
            PropertyDataType property_value =
                SymmTensor{values[0], values[1], values[2],
                           values[3], values[4], values[5]};
            return std::make_unique<Constant>(property_value);
        }
        case 9:
        {
            // Tensor
            PropertyDataType property_value =
                Tensor{values[0], values[1], values[2], values[3], values[4],
                       values[5], values[6], values[7], values[8]};
            return std::make_unique<Constant>(property_value);
        }

        default:
        {
            OGS_FATAL(
                "Creation of a constant property with %i components is not "
                "implemented.",
                values.size());
        }
    }

    PropertyDataType property_value;
    return std::make_unique<Constant>(property_value);
}
}  // namespace MaterialPropertyLib