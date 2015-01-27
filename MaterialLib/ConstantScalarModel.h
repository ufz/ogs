/*!
   \file  ConstantScalarModel.h
   \brief Declaration of class ConstantScalarModel.

   \author Wenqing Wang
   \date Jan 2015

   \copyright
    Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef CONSTANT_SCALAR_MODEL_H_
#define CONSTANT_SCALAR_MODEL_H_

#include<string>

#include"DensityType.h"

namespace MaterialLib
{

/// Constant scalar parameter.
class ConstantScalarModel
{
    public:
        ConstantScalarModel(const double value) : _value(value)
        {
        }

        /// Get model name.
        std::string getName() const
        {
            return "Constant";
        }

        unsigned getType() const
        {
            return 0;
        }

        /// Get model value
        double getValue() const
        {
            return _value;
        }
    private:
        double _value;
};

} // end namespace
#endif
