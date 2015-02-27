/*!
   \file  TensorParameter.h
   \brief Declaration of generic class for material parameters given in
          a tensor.

   \author Wenqing Wang
   \date Feb 2015

   \copyright
    Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef TENSOR_PARAMETER_H_
#define TENSOR_PARAMETER_H_

#include <memory>
#include <string>

#include <MaterialLib/ParameterBase.h>

namespace MaterialLib
{
/// Generic tensor parameter.
template<typename T_PARAMETER_TYPE, typename T_MAT_MODEL, typename T_MATRIX> class TensorParameter
    : public ParameterBase<T_PARAMETER_TYPE>
{
    public:

        template<typename... Args> TensorParameter(Args... args)
        {
            _material_model.reset(new T_MAT_MODEL(args...));
            this->setType(static_cast<T_PARAMETER_TYPE>(_material_model->getType()));
        }

        /// Get the model name.
        std::string getName() const
        {
            return _material_model->getName();
        }

        /// Get the parameter tensor in a matrix.
        template<typename... Args> T_MATRIX &getParameterMatrix(Args... args) const
        {
            return _material_model->getValue(args...);
        }
    private:
        std::unique_ptr<T_MAT_MODEL> _material_model;
};

} // end namespace
#endif
