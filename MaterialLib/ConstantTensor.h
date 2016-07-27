/*!
   \file  ConstantTensor.h
   \brief Declaration of constant tensor for material parameters.

   \author Wenqing Wang
   \date Feb 2015

   \copyright
    Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef CONSTANT_TENSOR_H_
#define CONSTANT_TENSOR_H_

namespace MaterialLib
{
/// Keep intrisinc permeability.
template<typename T_MATRIX> class ConstantTensor
{
    public:
        ConstantTensor(const T_MATRIX &perm) : _perm(perm)
        {
        }

        /// Get model name.
        std::string getName() const
        {
            return "Constant tensor parameter";
        }

        unsigned getType() const
        {
            return 1;
        }

        /// Get permeability tensor.
        T_MATRIX &getValue()
        {
            return _perm;
        }
    private:
        T_MATRIX _perm;
};

} // end namespace
#endif
