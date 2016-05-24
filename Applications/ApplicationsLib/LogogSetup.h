/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef APPLICATIONSLIB_LOGOGSETUP_H_
#define APPLICATIONSLIB_LOGOGSETUP_H_

#include <logog/include/logog.hpp>

#include <map>
#include <memory>
#include <string>

#include "BaseLib/LogogSimpleFormatter.h"

namespace ApplicationsLib
{

/// Initialization and shutting down of the logog library.
class LogogSetup final
{
public:

    LogogSetup()
    {
        LOGOG_INITIALIZE();
        logog_cout = std::unique_ptr<logog::Cout>(new logog::Cout);
        SetFormatter(std::unique_ptr<BaseLib::LogogSimpleFormatter>
            (new BaseLib::LogogSimpleFormatter));
    }

    ~LogogSetup()
    {
        // Objects have to be deleted before shutdown
        fmt.reset(nullptr);
        logog_cout.reset(nullptr);
        LOGOG_SHUTDOWN();
    }

    void SetFormatter(std::unique_ptr<logog::Formatter>&& formatter)
    {
        fmt = std::move(formatter);
        logog_cout->SetFormatter(*fmt);
    }

    void SetLevel(LOGOG_LEVEL_TYPE level)
    {
        logog::SetDefaultLevel(level);
    }

    void SetLevel(std::string const & level)
    {
        std::map<std::string, LOGOG_LEVEL_TYPE> foo =
        {
            { "none", LOGOG_LEVEL_NONE },
            { "emergency", LOGOG_LEVEL_EMERGENCY },
            { "alert", LOGOG_LEVEL_ALERT},
            { "critical", LOGOG_LEVEL_CRITICAL },
            { "error", LOGOG_LEVEL_ERROR },
            { "warn", LOGOG_LEVEL_WARN },
            { "info", LOGOG_LEVEL_INFO },
            { "debug", LOGOG_LEVEL_DEBUG },
            { "all", LOGOG_LEVEL_ALL }
        };


        //LOGOG_LEVEL_TYPE level_type;
        if(foo.find(level) != foo.end())
            SetLevel(foo[level]);
        else
        {
            ERR("%s is not a valid log level! Aborting.", level.c_str());
            std::abort();
        }
    }

private:
    std::unique_ptr<logog::Formatter> fmt;
    std::unique_ptr<logog::Cout> logog_cout;
};

}    // ApplicationsLib

#endif  // APPLICATIONSLIB_LOGOGSETUP_H_
