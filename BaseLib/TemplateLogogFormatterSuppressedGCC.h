/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file TemplateLogogFormatterSuppressedGCC.h
 *
 * Created on 2012-12-06 by Norihiro Watanabe
 */

#ifndef TEMPLATELOGOGFORMATTERSUPPRESSEDGCC_H_
#define TEMPLATELOGOGFORMATTERSUPPRESSEDGCC_H_

// ** INCLUDES **
#include "logog/include/logog.hpp"

namespace BaseLib {

/**
 * \brief TemplateLogogFormatterSuppressedGCC strips topics given as a template
 * parameter from logog::FormatterGCC.
 * See http://johnwbyrd.github.com/logog/customformatting.html for details.
 **/
template <int T_SUPPPRESS_TOPIC_FLAG>
class TemplateLogogFormatterSuppressedGCC : public logog::FormatterGCC
{

	virtual TOPIC_FLAGS GetTopicFlags( const logog::Topic &topic )
	{
	return ( logog::Formatter::GetTopicFlags( topic ) &
		~( T_SUPPPRESS_TOPIC_FLAG ));
	}

};

#endif // TEMPLATELOGOGFORMATTERSUPPRESSEDGCC_H_

} // namespace BaseLib
