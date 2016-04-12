/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_PROCESSOUTPUT_H
#define PROCESSLIB_PROCESSOUTPUT_H

namespace ProcessLib
{

template<typename GlobalVector>
struct SecondaryVariableFunctions
{
	using Fct = std::function<GlobalVector(
		GlobalVector const& x,
		AssemblerLib::LocalToGlobalIndexMap const& dof_table)>;

	Fct eval_field;
	Fct eval_residuals;
};

template<typename GlobalVector>
struct SecondaryVariable
{
	std::string const name;
	const unsigned n_components;
	SecondaryVariableFunctions<GlobalVector> fcts;
};

template <typename GlobalVector>
struct ProcessOutput
{
	std::vector<SecondaryVariable<GlobalVector>> secondary_variables;
	std::set<std::string> output_variables;

	bool output_residuals = false;
	//! Output global matrix/rhs after first iteration.
	bool output_global_matrix = false;
	bool output_iteration_results = false;
};

} // ProcessLib


#endif // PROCESSLIB_PROCESSOUTPUT_H
