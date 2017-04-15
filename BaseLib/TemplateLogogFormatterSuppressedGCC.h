/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-12-06
 * \brief  Definition of the TemplateLogogFormatterSuppressedGCC class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <string>

#include <logog/include/logog.hpp>
#ifdef USE_MPI
#include <mpi.h>
#endif

#include "StringTools.h"

namespace BaseLib {

/**
 * \brief TemplateLogogFormatterSuppressedGCC strips topics given as a template
 * parameter from logog::FormatterGCC.
 * See http://johnwbyrd.github.com/logog/customformatting.html for details.
 **/
template <int T_SUPPPRESS_TOPIC_FLAG>
class TemplateLogogFormatterSuppressedGCC : public logog::FormatterGCC
{
public:
#ifdef USE_MPI
    TemplateLogogFormatterSuppressedGCC(MPI_Comm mpi_comm = MPI_COMM_WORLD);
#endif

    TOPIC_FLAGS GetTopicFlags(const logog::Topic& topic) override
    {
    return ( logog::Formatter::GetTopicFlags( topic ) &
        ~( T_SUPPPRESS_TOPIC_FLAG ));
    }

    LOGOG_STRING& Format(const logog::Topic& topic,
                         const logog::Target& target) override;

private:
#ifdef USE_MPI
    std::string _str_mpi_rank;
#endif
};

} // namespace BaseLib

#include "TemplateLogogFormatterSuppressedGCC-impl.h"
