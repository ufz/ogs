/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file testrunner.cpp
 * Created on 2012-04-29 by Lars Bilke
 *
 */

// ** INCLUDES **
#include "gtest/gtest.h"
#include "logog/include/logog.hpp"
#include "logog/include/formatter.hpp"
#ifdef USE_LIS
#include "lis.h"
#endif

/**
 * new formatter for logog
 */
class FormatterCustom : public logog::FormatterGCC
{
    virtual TOPIC_FLAGS GetTopicFlags( const logog::Topic &topic )
    {
        return ( Formatter::GetTopicFlags( topic ) &
                 ~( TOPIC_LEVEL_FLAG | TOPIC_FILE_NAME_FLAG | TOPIC_LINE_NUMBER_FLAG ));
    }
};

/// Implementation of the googletest testrunner
int main(int argc, char* argv[])
{
    int ret = 0;
    LOGOG_INITIALIZE();
    logog::Cout out;
    FormatterCustom custom_format;
    out.SetFormatter(custom_format);

    try {
#ifdef USE_LIS
        lis_initialize(&argc, &argv);
#endif
        testing::InitGoogleTest ( &argc, argv );
        ret = RUN_ALL_TESTS();
    } catch (char* e) {
        ERR(e);
    } catch (std::exception& e) {
        ERR(e.what());
    } catch (...) {
        ERR("Unknown exception occurred!");
    }
#ifdef USE_LIS
    lis_finalize();
#endif
    LOGOG_SHUTDOWN();

    return ret;
}
