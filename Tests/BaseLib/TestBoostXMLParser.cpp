#include <gtest/gtest.h>
#include <logog/include/logog.hpp>

#include "BaseLib/ConfigTree.h"
#include "Tests/TestTools.h"

TEST(BaseLibBoostXML, T1)
{
    const char xml[] =
            "<double>5.6e-4</param>"
            ;
    auto const ptree = readXml(xml);
    BaseLib::ConfigTree ctree(ptree, "", BaseLib::ConfigTree::onerror,
                              BaseLib::ConfigTree::onwarning);
    auto const d = ctree.getConfigParameter<double>("double");
    INFO("double value: %g", d);
}
