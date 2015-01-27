/*!
   \file  ConstantDensityModel.h
   \brief Declaration of class ConstantDensityModel for constant density.

   \author Wenqing Wang
   \date Jan 2015

   \copyright
    Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef CONSTANT_DENSITY_MODEL_H_
#define CONSTANT_DENSITY_MODEL_H_

#include<string>

#include"DensityType.h"

namespace MaterialLib
{

/// Constant density model
class ConstantDensityModel
{
    public:
        ConstantDensityModel(const double density) : _density(density)
        {
        }

        /// Get desity model name.
        std::string getName() const
        {
            return "Constant";
        }

        DensityType getType() const
        {
            return DensityType::CONSTANT;
        }

        /// Get density value
        double getValue() const
        {
            return _density;
        }
    private:
        double _density;
};

} // end namespace
#endif

