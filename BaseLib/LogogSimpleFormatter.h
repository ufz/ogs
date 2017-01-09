/**
 * \file
 * \author Lars Bilke
 * \date   2012-09-13
 * \brief  Definition of the LogogSimpleFormatter class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef LOGOGSIMPLEFORMATTER_H
#define LOGOGSIMPLEFORMATTER_H

#include <logog/include/logog.hpp>

namespace BaseLib
{
/**
 * \brief LogogSimpleFormatter strips file name and line number from logog
 * output.
 * See http://johnwbyrd.github.com/logog/customformatting.html for details.
 **/
class LogogSimpleFormatter : public logog::FormatterMSVC
{
    virtual TOPIC_FLAGS GetTopicFlags(const logog::Topic& topic)
    {
        return (logog::Formatter::GetTopicFlags(topic) &
                ~(TOPIC_FILE_NAME_FLAG | TOPIC_LINE_NUMBER_FLAG));
    }
};

}  // namespace BaseLib

#endif  // LOGOGSIMPLEFORMATTER_H
