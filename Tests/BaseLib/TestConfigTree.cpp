/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>
#include <logog/include/logog.hpp>

#include <boost/property_tree/xml_parser.hpp>
#include <sstream>

#include "BaseLib/ConfigTreeNew.h"

#define DO_EXPECT(cbs, error, warning) do  { \
        if (error)   EXPECT_TRUE((cbs).get_error());   else EXPECT_FALSE((cbs).get_error()); \
        if (warning) EXPECT_TRUE((cbs).get_warning()); else EXPECT_FALSE((cbs).get_warning()); \
        (cbs).reset(); \
    } while(false)

#define RUN_SAFE(expr) do { \
        try { expr; } catch(Exc) {} \
    } while (false)

// Exception thrown by the error callback of the class below
class Exc {};

// class that provides callback functions used with ConfigTreeNew
class Callbacks
{
public:
    BaseLib::ConfigTreeNew::Callback
    get_error_cb() {
        return [this](std::string const& path, std::string const& message)
        {
            (void) path; (void) message;
            DBUG("error <%s> : %s", path.c_str(), message.c_str());
            _error = true;
            throw Exc(); // throw in order to stop normal execution
        };
    }

    BaseLib::ConfigTreeNew::Callback
    get_warning_cb() {
        return [this](std::string const& path, std::string const& message)
        {
            (void) path; (void) message;
            DBUG("warning <%s> : %s", path.c_str(), message.c_str());
            _warning = true;
        };
    }

    bool get_error()   const { return _error; }
    bool get_warning() const { return _warning; }
    void reset() { _error = false; _warning = false; }

private:
    bool _error = false;
    bool _warning = false;
};


TEST(BaseLibConfigTree, ConfigTreeEmpty)
{
    boost::property_tree::ptree ptree;
    Callbacks cbs;

    {
        BaseLib::ConfigTreeNew conf(ptree, cbs.get_error_cb(), cbs.get_warning_cb());
        (void) conf;
    } // ConfigTree destroyed here

    DO_EXPECT(cbs, false, false);
}


TEST(BaseLibConfigTree, ConfigTreeGet)
{
    const char xml[] =
            "<double>5.6e-4</double>"
            "<bool>true</bool>"
            "<int>5</int>"
            "<sub>"
            "  <float>6.1</float>"
            "  <float2>0.1</float2>"
            "  <ignored/>"
            "  <ignored2/>"
            "  <ignored2/>"
            "</sub>"
            "<x>Y</x>";

    boost::property_tree::ptree ptree;
    std::istringstream xml_str(xml);
    read_xml(xml_str, ptree);

    Callbacks cbs;
    {
        BaseLib::ConfigTreeNew conf(ptree, cbs.get_error_cb(), cbs.get_warning_cb());

        EXPECT_EQ(5.6e-4, conf.getConfParam<double>("double")); // read certain types
        DO_EXPECT(cbs, false, false);
        EXPECT_TRUE(conf.getConfParam<bool>("bool"));
        DO_EXPECT(cbs, false, false);
        EXPECT_EQ(5, conf.getConfParam<int>("int"));
        DO_EXPECT(cbs, false, false);

        EXPECT_EQ(8, conf.getConfParam<int>("intx", 8)); // reading with default value
        DO_EXPECT(cbs, false, false);

        {
            auto sub = conf.getConfSubtree("sub");
            DO_EXPECT(cbs, false, false);

            EXPECT_EQ(6.1f, sub.getConfParam<float>("float"));
            DO_EXPECT(cbs, false, false);

            if (auto f2 = sub.getConfParamOptional<float>("float2")) { // read optional value
                EXPECT_EQ(0.1f, *f2);
                DO_EXPECT(cbs, false, false);
            }

            auto f3 = sub.getConfParamOptional<float>("float3"); // optional value not existent
            ASSERT_FALSE(f3);
            DO_EXPECT(cbs, false, false);

            sub.ignoreConfParam("ignored");
            DO_EXPECT(cbs, false, false);
            sub.ignoreConfParamAll("ignored2");
            DO_EXPECT(cbs, false, false);
            sub.ignoreConfParamAll("ignored4"); // I can ignore nonexistent stuff
            DO_EXPECT(cbs, false, false);

            // I can not ignore stuff that I already read
            // this also makes sure that the subtree inherits the callbacks properly
            RUN_SAFE(sub.ignoreConfParam("float"));
            DO_EXPECT(cbs, true, false);
        }
        for (int i : {0, 1, 2}) {
            (void) i;
            EXPECT_EQ("Y", conf.peekConfParam<std::string>("x"));
            DO_EXPECT(cbs, false, false);
        }
        conf.checkConfParam<std::string>("x", "Y");
        DO_EXPECT(cbs, false, false);
    } // ConfigTree destroyed here
    DO_EXPECT(cbs, false, false);
}


TEST(BaseLibConfigTree, ConfigTreeIncompleteParse)
{
    const char xml[] =
            "<double>5.6</double>"
            "<bool>true</bool>"
            ;

    boost::property_tree::ptree ptree;
    std::istringstream xml_str(xml);
    read_xml(xml_str, ptree);

    Callbacks cbs;
    {
        BaseLib::ConfigTreeNew conf(ptree, cbs.get_error_cb(), cbs.get_warning_cb());

        EXPECT_EQ(5.6, conf.getConfParam<double>("double")); // read certain types
        DO_EXPECT(cbs, false, false);
    } // ConfigTree destroyed here
    DO_EXPECT(cbs, false, true); // expect warning because I didn't read everything
}


TEST(BaseLibConfigTree, ConfigTreeCheckRange)
{
    const char xml[] =
            "<val><int>0</int></val>"
            "<val><int>1</int></val>"
            "<val><int>2</int></val>"
            "<int>0</int>"
            "<int>1</int>"
            "<int>2</int>";

    boost::property_tree::ptree ptree;
    std::istringstream xml_str(xml);
    read_xml(xml_str, ptree);

    Callbacks cbs;
    {
        BaseLib::ConfigTreeNew conf(ptree, cbs.get_error_cb(), cbs.get_warning_cb());

        {
            auto list = conf.getConfSubtreeList("val");
            DO_EXPECT(cbs, false, false);
            EXPECT_EQ(3, std::distance(list.begin(), list.end()));
            DO_EXPECT(cbs, false, false);
            EXPECT_EQ(3, std::distance(list.begin(), list.end()));
            DO_EXPECT(cbs, false, false);
        }

        {
            auto list = conf.getConfParamList<int>("int");
            DO_EXPECT(cbs, false, false);
            EXPECT_EQ(3, std::distance(list.begin(), list.end()));
            DO_EXPECT(cbs, false, false);
            EXPECT_EQ(3, std::distance(list.begin(), list.end()));
            DO_EXPECT(cbs, false, false);
        }

    } // ConfigTree destroyed here

    // there will be warnings because I don't process the list entries
    DO_EXPECT(cbs, false, true);
}


TEST(BaseLibConfigTree, ConfigTreeGetSubtreeList)
{
    const char xml[] =
            "<val><int>0</int></val>"
            "<val><int>1</int></val>"
            "<val><int>2</int></val>";

    boost::property_tree::ptree ptree;
    std::istringstream xml_str(xml);
    read_xml(xml_str, ptree);

    Callbacks cbs;
    {
        BaseLib::ConfigTreeNew conf(ptree, cbs.get_error_cb(), cbs.get_warning_cb());

        int i = 0;
        for (auto ct : conf.getConfSubtreeList("val"))
        {
            EXPECT_EQ(i, ct.getConfParam<int>("int"));
            DO_EXPECT(cbs, false, false);
            ++i;
        }
    } // ConfigTree destroyed here
    DO_EXPECT(cbs, false, false);
}


TEST(BaseLibConfigTree, ConfigTreeGetValueList)
{
    const char xml[] =
            "<int>0</int>"
            "<int>1</int>"
            "<int>2</int>";

    boost::property_tree::ptree ptree;
    std::istringstream xml_str(xml);
    read_xml(xml_str, ptree);

    Callbacks cbs;
    {
        BaseLib::ConfigTreeNew conf(ptree, cbs.get_error_cb(), cbs.get_warning_cb());

        int n = 0;
        for (auto i : conf.getConfParamList<int>("int"))
        {
            EXPECT_EQ(n, i);
            DO_EXPECT(cbs, false, false);
            ++n;
        }
    } // ConfigTree destroyed here
    DO_EXPECT(cbs, false, false);
}


TEST(BaseLibConfigTree, ConfigTreeNoConversion)
{
    const char xml[] =
            "<int>5.6</int>"         // not convertible to int
            "<double>5.6tz</double>" // not convertible to double
            "<non_double>0.1x</non_double>" // not either convertible to double
            "<bool>true</bool>"
            "<ign/>"
            "<ign2/><ign2/><ign2/>"
            ;

    boost::property_tree::ptree ptree;
    std::istringstream xml_str(xml);
    read_xml(xml_str, ptree);

    Callbacks cbs;
    {
        BaseLib::ConfigTreeNew conf(ptree, cbs.get_error_cb(), cbs.get_warning_cb());

        RUN_SAFE(conf.getConfParam<int>("int"));
        DO_EXPECT(cbs, true, false);
        RUN_SAFE(conf.ignoreConfParam("int")); // after failure I also cannot ignore something
        DO_EXPECT(cbs, true, false);

        RUN_SAFE(conf.getConfParam<double>("double"));
        DO_EXPECT(cbs, true, false);

        // peek value existent but not convertible
        RUN_SAFE(conf.peekConfParam<double>("non_double"));
        DO_EXPECT(cbs, true, false);

        // optional value existent but not convertible
        RUN_SAFE(
            auto d = conf.getConfParamOptional<double>("non_double");
            ASSERT_FALSE(d);
        );
        DO_EXPECT(cbs, true, false);

        // assert that I can only ignore something once
        RUN_SAFE(conf.ignoreConfParam("ign"));
        DO_EXPECT(cbs, false, false);
        RUN_SAFE(conf.ignoreConfParam("ign"));
        DO_EXPECT(cbs, true, false);
        RUN_SAFE(conf.ignoreConfParamAll("ign2"));
        DO_EXPECT(cbs, false, false);
        RUN_SAFE(conf.ignoreConfParamAll("ign2"));
        DO_EXPECT(cbs, true, false);

        // assert that I cannot read a parameter twice
        RUN_SAFE(conf.getConfParam<bool>("bool"));
        DO_EXPECT(cbs, false, false);
        RUN_SAFE(conf.getConfParam<bool>("bool"));
        DO_EXPECT(cbs, true, false);

    } // ConfigTree destroyed here

    // There will bewarnings because I don't succeed in reading every setting,
    // and furthermore I read some setting too often.
    DO_EXPECT(cbs, false, false);
}


TEST(BaseLibConfigTree, BadKeynames)
{
    const char xml[] = "";

    boost::property_tree::ptree ptree;
    std::istringstream xml_str(xml);
    read_xml(xml_str, ptree);

    Callbacks cbs;
    {
        BaseLib::ConfigTreeNew conf(ptree, cbs.get_error_cb(), cbs.get_warning_cb());

        for (std::string tag : { "<", "Z", ".", "$", "0", "", "/" })
        {
            RUN_SAFE(conf.getConfParam<int>(tag));
            DO_EXPECT(cbs, true, false);
            RUN_SAFE(conf.getConfParam<int>(tag, 500));
            DO_EXPECT(cbs, true, false);
            RUN_SAFE(conf.getConfParamOptional<int>(tag));
            DO_EXPECT(cbs, true, false);
            RUN_SAFE(conf.getConfParamList<int>(tag));

            DO_EXPECT(cbs, true, false);
            RUN_SAFE(conf.peekConfParam<int>(tag));
            DO_EXPECT(cbs, true, false);
            RUN_SAFE(conf.checkConfParam<int>(tag, 500));

            DO_EXPECT(cbs, true, false);
            RUN_SAFE(conf.getConfSubtree(tag));
            DO_EXPECT(cbs, true, false);
            RUN_SAFE(conf.getConfSubtreeOptional(tag));
            DO_EXPECT(cbs, true, false);
            RUN_SAFE(conf.getConfSubtreeList(tag));
            DO_EXPECT(cbs, true, false);
        }

    } // ConfigTree destroyed here

    DO_EXPECT(cbs, false, false);
}

// String literals are somewhat special for template classes
TEST(BaseLibConfigTree, ConfigTreeStringLiterals)
{
    const char xml[] =
            "<s>test</s>"
            "<t>Test</t>";

    boost::property_tree::ptree ptree;
    std::istringstream xml_str(xml);
    read_xml(xml_str, ptree);

    Callbacks cbs;
    {
        BaseLib::ConfigTreeNew conf(ptree, cbs.get_error_cb(), cbs.get_warning_cb());

        EXPECT_EQ("test", conf.getConfParam<std::string>("s", "XX"));
        DO_EXPECT(cbs, false, false);

        EXPECT_EQ("XX",   conf.getConfParam<std::string>("n", "XX"));
        DO_EXPECT(cbs, false, false);

        conf.checkConfParam("t", "Test");
        DO_EXPECT(cbs, false, false);
    } // ConfigTree destroyed here
    DO_EXPECT(cbs, false, false);
}

// String literals are somewhat special for template classes
TEST(BaseLibConfigTree, ConfigTreeMove)
{
    const char xml[] =
            "<s>test</s>"
            "<t>Test</t>";

    boost::property_tree::ptree ptree;
    std::istringstream xml_str(xml);
    read_xml(xml_str, ptree);

    Callbacks cbs;
    {
        BaseLib::ConfigTreeNew conf(ptree, cbs.get_error_cb(), cbs.get_warning_cb());

        EXPECT_EQ("test", conf.getConfParam<std::string>("s", "XX"));
        DO_EXPECT(cbs, false, false);

        BaseLib::ConfigTreeNew conf2(std::move(conf));

        EXPECT_EQ("XX",   conf2.getConfParam<std::string>("n", "XX"));
        DO_EXPECT(cbs, false, false);

        conf2.checkConfParam("t", "Test");
        DO_EXPECT(cbs, false, false);
    } // ConfigTree destroyed here
    DO_EXPECT(cbs, false, false);
}

