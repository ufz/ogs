/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-06-25
 * \brief  Definition of the LisOption class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef LIS_OPTION_H_
#define LIS_OPTION_H_

#include <string>
#include <map>

#include "BaseLib/ConfigTree.h"

namespace MathLib
{

/**
 * \brief Option for Lis solver
 */
struct LisOption
{
    using Key = std::string;
    using Value = std::string;

    std::map<Key, Value> settings;

    LisOption()
    {
        // by default do not zero out the solution vector
        settings["-initxzeros"] = "0";
    }

    /**
     * Adds options from a \c ConfigTree.
     *
     * Options are stored as strings as obtained from the ConfigTree.
     * In particular it will not be checked by this method if options are valid.
     *
     * It is only guaranteed that each given option will be passed only once to Lis.
     *
     * The \c ConfigTree section must have the layout as shown in the
     * following example:
     *
     * \code{.xml}
     * <linear_solver>
     *   <solver_type>         gmres </solver_type>
     *   <precon_type>           ilu </precon_type>
     *   <error_tolerance>   1.0e-10 </error_tolerance>
     *   <max_iteration_step>    200 </max_iteration_step>
     *
     *   <lis_option><option> -ilu_fill </option><value>   3 </value></lis_option>
     *   <lis_option><option>    -print </option><value> mem </value></lis_option>
     * </linear_solver>
     * \endcode
     *
     * The first four tags are special in the sense that those options
     * will take precedence even if there is an equivalent Lis option specified.
     *
     * Any possible Lis option can be supplied by a line like
     *
     * \code{.xml}
     *   <lis_option><option> -ilu_fill </option><value>   3 </value></lis_option>
     * \endcode
     *
     * For possible values, please refer to the Lis User Guide:
     * http://www.ssisc.org/lis/lis-ug-en.pdf
     */
    void addOptions(BaseLib::ConfigTree const& options);

    void printInfo() const;


    /// Matrix type
    enum class MatrixType : int
    {
        CRS = 1,
        CCS = 2,
        MSR = 3,
        DIA = 4,
        ELL = 5,
        JDS = 6,
        BSR = 7,
        BSC = 8,
        VBR = 9,
        COO = 10,
        DNS = 11
    };
};

}
#endif //LIS_OPTION_H_
