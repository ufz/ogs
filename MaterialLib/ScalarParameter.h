/*!
   \file  ScalarParameter.h
   \brief Declaration of generic class for scalar material parameter.

   \author Wenqing Wang
   \date Jan 2015

   \copyright
    Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
 */
#ifndef SCALAR_PARAMETER_H_
#define SCALAR_PARAMETER_H_

#include <memory>
#include <string>

#include "MaterialLib/ParameterBase.h"

namespace MaterialLib
{
/// Class for scalar material parameter
template<typename T_PARAMETER_TYPE,typename T_MAT_MODEL> class ScalarParameter
: public ParameterBase<T_PARAMETER_TYPE>
{
public:

    template<typename... Args> ScalarParameter(Args... args)
    {
        _material_model.reset(new T_MAT_MODEL(args...));
        this->setType(static_cast<T_PARAMETER_TYPE> (_material_model->getType()));
    }

    /// Get material model name.
    std::string getName() const
    {
        return _material_model->getName();
    }

    /// Get the parameter value.
    template<typename... Args> double getValue(Args... args) const
    {
        return _material_model->getValue(args...);
    }

    /// Get the partial differential of the parameter with respect to variables.
    template<typename... Args> double getdValue(Args... args) const
    {
        return _material_model->getdValue(args...);
    }

private:
    std::unique_ptr<T_MAT_MODEL> _material_model;
};

} // end namespace
#endif
