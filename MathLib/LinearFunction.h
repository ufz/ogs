/**
 * \author Norihiro Watanabe
 * \date   2013-08-13
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef LINEARFUNCTION_H_
#define LINEARFUNCTION_H_

#include <array>
#include <algorithm>
#include <cassert>

namespace MathLib
{

/**
 * Linear function
 * \f[
 * 	f(x_1,...,x_k)=a_0+a_1*x_1+...+a_k*x_k
 * \f]
 *
 * \tparam T_TYPE  value type
 * \tparam N_VARS  the number of variables (k)
 */
template <typename T_TYPE, unsigned N_VARS>
class LinearFunction
{
public:
	/**
	 * Constructor
	 * \param coefficients  an array of coefficients of a linear function.
	 * The size of the coefficient array should equal to the number of variables + 1
	 */
	explicit LinearFunction(const T_TYPE* coefficients)
	: _coefficients(copy(coefficients))
	{
		static_assert(N_VARS>0, "Template parameter N_VARS in LinearFunction should be nonzero");
	}

	/**
	 * evaluate the function
	 * \param x  an array of variables. the size of the array should equal to the number of variables + 1
	 */
	T_TYPE operator()(T_TYPE const * const x) const
	{
		T_TYPE v = _coefficients[0];
		for (unsigned i=0; i<N_VARS; i++)
			v += _coefficients[i+1]*x[i];
		return v;
	}

private:
	/// copy a given raw array to std::array
	static std::array<T_TYPE, N_VARS+1> copy(const T_TYPE* coefficients)
	{
		std::array<T_TYPE, N_VARS+1> temp;
		std::copy(coefficients, coefficients+ temp.size(), temp.data());
		return temp;
	}

	/// Coefficients of a linear function
	const std::array<T_TYPE, N_VARS+1> _coefficients;
};

} // MathLib

#endif /* LINEARFUNCTION_H_ */
