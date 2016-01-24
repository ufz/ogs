/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
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

// make useful line numbers appear in the output of gtest
#define EXPECT_ERR_WARN(cbs, error, warning) do  { \
        if (error)   EXPECT_TRUE((cbs).get_error());   else EXPECT_FALSE((cbs).get_error()); \
        if (warning) EXPECT_TRUE((cbs).get_warning()); else EXPECT_FALSE((cbs).get_warning()); \
        (cbs).reset(); \
    } while(false)

// catch exception thrown by the error callback
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
    get_error_cb()
    {
        return [this](std::string const& path, std::string const& message)
        {
            (void) path; (void) message;
            DBUG("error <%s> : %s", path.c_str(), message.c_str());
            _error = true;
            throw Exc(); // throw in order to stop normal execution
        };
    }

    BaseLib::ConfigTreeNew::Callback
    get_warning_cb()
    {
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


boost::property_tree::ptree
readXml(const char xml[])
{
    boost::property_tree::ptree ptree;
    std::istringstream xml_str(xml);
    read_xml(xml_str, ptree,
             boost::property_tree::xml_parser::no_comments |
             boost::property_tree::xml_parser::trim_whitespace);
    return ptree;
}


TEST(BaseLibConfigTree, Empty)
{
    boost::property_tree::ptree ptree;
    Callbacks cbs;

    {
        BaseLib::ConfigTreeNew conf(ptree, cbs.get_error_cb(), cbs.get_warning_cb());
        (void) conf;
    } // ConfigTree destroyed here

    EXPECT_ERR_WARN(cbs, false, false);
}


TEST(BaseLibConfigTree, Get)
{
    const char xml[] =
            "<double>5.6e-4</double>"
            "<bool>true</bool>"
            "<int>5</int>"
            "<sub>"
            "  <float>6.1</float>"
            "  <float2>0.1</float2>"
            "  <bool1>false</bool1>"
            "  <bool2>false</bool2>"
            "  <bool3/>"
            "  <ignored/>"
            "  <ignored2/>"
            "  <ignored2/>"
            "</sub>"
            "<x>Y</x>"
            "<z attr=\"0.5\">32.0</z>";
    auto const ptree = readXml(xml);

    Callbacks cbs;
    {
        BaseLib::ConfigTreeNew conf(ptree, cbs.get_error_cb(), cbs.get_warning_cb());

        EXPECT_EQ(5.6e-4, conf.getConfParam<double>("double")); // read certain types
        EXPECT_ERR_WARN(cbs, false, false);
        EXPECT_TRUE(conf.getConfParam<bool>("bool"));
        EXPECT_ERR_WARN(cbs, false, false);
        EXPECT_EQ(5, conf.getConfParam<int>("int"));
        EXPECT_ERR_WARN(cbs, false, false);

        EXPECT_EQ(8, conf.getConfParam<int>("intx", 8)); // reading with default value
        EXPECT_ERR_WARN(cbs, false, false);


        // Testing subtree
        {
            auto sub = conf.getConfSubtree("sub");
            EXPECT_ERR_WARN(cbs, false, false);

            EXPECT_EQ(6.1f, sub.getConfParam<float>("float"));
            EXPECT_ERR_WARN(cbs, false, false);

            if (auto f2 = sub.getConfParamOptional<float>("float2")) { // read optional value
                EXPECT_EQ(0.1f, *f2);
            }
            EXPECT_ERR_WARN(cbs, false, false);

            auto f3 = sub.getConfParamOptional<float>("float3"); // optional value not existent
            ASSERT_FALSE(f3);
            EXPECT_ERR_WARN(cbs, false, false);


            // Testing the getConfParam...() (non-template) / getValue() combination

            auto bool1 = sub.getConfParam("bool1");
            EXPECT_ERR_WARN(cbs, false, false);
            EXPECT_FALSE(bool1.getValue<bool>());
            EXPECT_ERR_WARN(cbs, false, false);
            RUN_SAFE(bool1.getValue<bool>()); // getting data twice
            EXPECT_ERR_WARN(cbs, true, false);

            if (auto bool2 = sub.getConfParamOptional("bool2")) {
                EXPECT_ERR_WARN(cbs, false, false);
                EXPECT_FALSE(bool2->getValue<bool>());
            }
            EXPECT_ERR_WARN(cbs, false, false);

            if (auto bool3 = sub.getConfParamOptional("bool3")) {
                EXPECT_ERR_WARN(cbs, false, false);
                RUN_SAFE(bool3->getValue<bool>());
                EXPECT_ERR_WARN(cbs, true, false); // error because of no data
            }
            EXPECT_ERR_WARN(cbs, false, false);

            EXPECT_FALSE(sub.getConfParamOptional("bool4")); // optional value not existent
            EXPECT_ERR_WARN(cbs, false, false);


            // Testing ignore

            sub.ignoreConfParam("ignored");
            EXPECT_ERR_WARN(cbs, false, false);
            sub.ignoreConfParamAll("ignored2");
            EXPECT_ERR_WARN(cbs, false, false);
            sub.ignoreConfParamAll("ignored4"); // I can ignore nonexistent stuff
            EXPECT_ERR_WARN(cbs, false, false);

            // I can not ignore stuff that I already read
            // this also makes sure that the subtree inherits the callbacks properly
            RUN_SAFE(sub.ignoreConfParam("float"));
            EXPECT_ERR_WARN(cbs, true, false);
        }
        for (int i : {0, 1, 2}) {
            (void) i;
            EXPECT_EQ("Y", conf.peekConfParam<std::string>("x"));
            EXPECT_ERR_WARN(cbs, false, false);
        }
        conf.checkConfParam<std::string>("x", "Y");
        EXPECT_ERR_WARN(cbs, false, false);


        // Testing attributes
        {
            auto z = conf.getConfParam("z");
            EXPECT_ERR_WARN(cbs, false, false);
            EXPECT_EQ(0.5, z.getConfAttribute<double>("attr"));
            EXPECT_ERR_WARN(cbs, false, false);
            RUN_SAFE(z.getConfAttribute<double>("attr")); // getting attr twice
            EXPECT_ERR_WARN(cbs, true, false);
            RUN_SAFE(z.getConfAttribute<double>("not_an_attr")); // nonexistent attribute
            EXPECT_ERR_WARN(cbs, true, false);
            EXPECT_EQ(32.0, z.getValue<double>());
            EXPECT_ERR_WARN(cbs, false, false);
        }
        EXPECT_ERR_WARN(cbs, false, false);
    } // ConfigTree destroyed here
    EXPECT_ERR_WARN(cbs, false, false);
}


TEST(BaseLibConfigTree, IncompleteParse)
{
    const char xml[] =
            "<double>5.6</double>"
            "<not_read>true</not_read>"
            "<tag>this data won't be read</tag>"
            "<pt x=\"0.5\">1</pt>"
            "<pt2 x=\"0.5\" y=\"1.0\" z=\"2.0\" />"
            ;
    auto const ptree = readXml(xml);

    Callbacks cbs;
    {
        BaseLib::ConfigTreeNew conf(ptree, cbs.get_error_cb(), cbs.get_warning_cb());

        EXPECT_EQ(5.6, conf.getConfParam<double>("double"));
        EXPECT_ERR_WARN(cbs, false, false);

        conf.getConfSubtree("tag");
        EXPECT_ERR_WARN(cbs, false, true); // data of <tag> has not been read

        EXPECT_EQ(1, conf.getConfParam<int>("pt"));
        EXPECT_ERR_WARN(cbs, false, true); // attribute "x" has not been read

        {
            auto pt2 = conf.getConfParam("pt2");
            EXPECT_EQ(0.5, pt2.getConfAttribute<double>("x"));
            EXPECT_ERR_WARN(cbs, false, false);
            EXPECT_EQ(1.0, pt2.getConfAttribute<double>("y"));
            EXPECT_ERR_WARN(cbs, false, false);

            BaseLib::checkAndInvalidate(pt2);
            EXPECT_ERR_WARN(cbs, false, true); // attribute "z" not read
        }
        EXPECT_ERR_WARN(cbs, false, false);

    } // ConfigTree destroyed here
    EXPECT_ERR_WARN(cbs, false, true); // expect warning because I didn't read everything
}


TEST(BaseLibConfigTree, CheckRange)
{
    const char xml[] =
            "<val><int>0</int></val>"
            "<val><int>1</int></val>"
            "<val><int>2</int></val>"
            "<int>0</int>"
            "<int>1</int>"
            "<int>2</int>";
    auto const ptree = readXml(xml);

    Callbacks cbs;
    {
        BaseLib::ConfigTreeNew conf(ptree, cbs.get_error_cb(), cbs.get_warning_cb());

        {
            // check that std::distance can be computed twice in a row
            auto list = conf.getConfSubtreeList("val");
            EXPECT_ERR_WARN(cbs, false, false);
            EXPECT_EQ(3, std::distance(list.begin(), list.end()));
            EXPECT_ERR_WARN(cbs, false, false);
            EXPECT_EQ(3, std::distance(list.begin(), list.end()));
            EXPECT_ERR_WARN(cbs, false, false);
        }

        {
            // check that std::distance can be computed twice in a row
            auto list = conf.getConfParamList<int>("int");
            EXPECT_ERR_WARN(cbs, false, false);
            EXPECT_EQ(3, std::distance(list.begin(), list.end()));
            EXPECT_ERR_WARN(cbs, false, false);
            EXPECT_EQ(3, std::distance(list.begin(), list.end()));
            EXPECT_ERR_WARN(cbs, false, false);
        }

    } // ConfigTree destroyed here

    // there will be warnings because I don't process the list entries
    EXPECT_ERR_WARN(cbs, false, true);
}


TEST(BaseLibConfigTree, GetSubtreeList)
{
    const char xml[] =
            "<val><int>0</int></val>"
            "<val><int>1</int></val>"
            "<val><int>2</int></val>";
    auto const ptree = readXml(xml);

    Callbacks cbs;
    {
        BaseLib::ConfigTreeNew conf(ptree, cbs.get_error_cb(), cbs.get_warning_cb());

        int i = 0;
        for (auto ct : conf.getConfSubtreeList("val"))
        {
            EXPECT_EQ(i, ct.getConfParam<int>("int"));
            EXPECT_ERR_WARN(cbs, false, false);
            ++i;
        }
    } // ConfigTree destroyed here
    EXPECT_ERR_WARN(cbs, false, false);
}


TEST(BaseLibConfigTree, GetValueList)
{
    const char xml[] =
            "<int>0</int>"
            "<int>1</int>"
            "<int>2</int>";
    auto const ptree = readXml(xml);

    Callbacks cbs;
    {
        BaseLib::ConfigTreeNew conf(ptree, cbs.get_error_cb(), cbs.get_warning_cb());

        int n = 0;
        for (auto i : conf.getConfParamList<int>("int"))
        {
            EXPECT_EQ(n, i);
            EXPECT_ERR_WARN(cbs, false, false);
            ++n;
        }
    } // ConfigTree destroyed here
    EXPECT_ERR_WARN(cbs, false, false);
}


TEST(BaseLibConfigTree, NoConversion)
{
    const char xml[] =
            "<int>5.6</int>"         // not convertible to int
            "<double>5.6tz</double>" // not convertible to double
            "<non_double>0.1x</non_double>" // not either convertible to double
            "<bool>true</bool>"
            "<ign/>"
            "<ign2/><ign2/><ign2/>"
            ;
    auto const ptree = readXml(xml);

    Callbacks cbs;
    {
        BaseLib::ConfigTreeNew conf(ptree, cbs.get_error_cb(), cbs.get_warning_cb());

        RUN_SAFE(conf.getConfParam<int>("int"));
        EXPECT_ERR_WARN(cbs, true, false);
        RUN_SAFE(conf.ignoreConfParam("int")); // after failure I also cannot ignore something
        EXPECT_ERR_WARN(cbs, true, false);

        RUN_SAFE(conf.getConfParam<double>("double"));
        EXPECT_ERR_WARN(cbs, true, false);

        // peek value existent but not convertible
        RUN_SAFE(conf.peekConfParam<double>("non_double"));
        EXPECT_ERR_WARN(cbs, true, false);

        // optional value existent but not convertible
        RUN_SAFE(
            auto d = conf.getConfParamOptional<double>("non_double");
            ASSERT_FALSE(d);
        );
        EXPECT_ERR_WARN(cbs, true, false);

        // assert that I can only ignore something once
        RUN_SAFE(conf.ignoreConfParam("ign"));
        EXPECT_ERR_WARN(cbs, false, false);
        RUN_SAFE(conf.ignoreConfParam("ign"));
        EXPECT_ERR_WARN(cbs, true, false);
        RUN_SAFE(conf.ignoreConfParamAll("ign2"));
        EXPECT_ERR_WARN(cbs, false, false);
        RUN_SAFE(conf.ignoreConfParamAll("ign2"));
        EXPECT_ERR_WARN(cbs, true, false);

        // assert that I cannot read a parameter twice
        RUN_SAFE(conf.getConfParam<bool>("bool"));
        EXPECT_ERR_WARN(cbs, false, false);
        RUN_SAFE(conf.getConfParam<bool>("bool"));
        EXPECT_ERR_WARN(cbs, true, false);

    } // ConfigTree destroyed here

    // There will bewarnings because I don't succeed in reading every setting,
    // and furthermore I read some setting too often.
    EXPECT_ERR_WARN(cbs, false, false);
}


TEST(BaseLibConfigTree, BadKeynames)
{
    const char xml[] = "";
    auto const ptree = readXml(xml);

    Callbacks cbs;
    {
        BaseLib::ConfigTreeNew conf(ptree, cbs.get_error_cb(), cbs.get_warning_cb());

        for (auto tag : { "<", "Z", ".", "$", "0", "", "/", "_", "a__" })
        {
            RUN_SAFE(conf.getConfParam<int>(tag));
            EXPECT_ERR_WARN(cbs, true, false);
            RUN_SAFE(conf.getConfParam<int>(tag, 500));
            EXPECT_ERR_WARN(cbs, true, false);
            RUN_SAFE(conf.getConfParamOptional<int>(tag));
            EXPECT_ERR_WARN(cbs, true, false);
            RUN_SAFE(conf.getConfParamList<int>(tag));
            EXPECT_ERR_WARN(cbs, true, false);

            RUN_SAFE(conf.getConfParam(tag));
            EXPECT_ERR_WARN(cbs, true, false);
            RUN_SAFE(conf.getConfParamOptional(tag));
            EXPECT_ERR_WARN(cbs, true, false);

            RUN_SAFE(conf.peekConfParam<int>(tag));
            EXPECT_ERR_WARN(cbs, true, false);
            RUN_SAFE(conf.checkConfParam<int>(tag, 500));
            EXPECT_ERR_WARN(cbs, true, false);

            RUN_SAFE(conf.getConfSubtree(tag));
            EXPECT_ERR_WARN(cbs, true, false);
            RUN_SAFE(conf.getConfSubtreeOptional(tag));
            EXPECT_ERR_WARN(cbs, true, false);
            RUN_SAFE(conf.getConfSubtreeList(tag));
            EXPECT_ERR_WARN(cbs, true, false);

            RUN_SAFE(conf.getConfAttribute<int>(tag));
            EXPECT_ERR_WARN(cbs, true, false);
        }

    } // ConfigTree destroyed here

    EXPECT_ERR_WARN(cbs, false, false);
}

// String literals are somewhat special for template classes
TEST(BaseLibConfigTree, StringLiterals)
{
    const char xml[] =
            "<s>test</s>"
            "<t>Test</t>";
    auto const ptree = readXml(xml);

    Callbacks cbs;
    {
        BaseLib::ConfigTreeNew conf(ptree, cbs.get_error_cb(), cbs.get_warning_cb());

        EXPECT_EQ("test", conf.getConfParam<std::string>("s", "XX"));
        EXPECT_ERR_WARN(cbs, false, false);

        // <n> not present in the XML, so return the default value
        EXPECT_EQ("XX",   conf.getConfParam<std::string>("n", "XX"));
        EXPECT_ERR_WARN(cbs, false, false);

        conf.checkConfParam("t", "Test");
        EXPECT_ERR_WARN(cbs, false, false);
    } // ConfigTree destroyed here
    EXPECT_ERR_WARN(cbs, false, false);
}

// String literals are somewhat special for template classes
TEST(BaseLibConfigTree, MoveConstruct)
{
    const char xml[] =
            "<s>test</s>"
            "<t>Test</t>"
            "<u>data</u>";
    auto const ptree = readXml(xml);

    Callbacks cbs;
    {
        BaseLib::ConfigTreeNew conf(ptree, cbs.get_error_cb(), cbs.get_warning_cb());

        EXPECT_EQ("test", conf.getConfParam<std::string>("s", "XX"));
        EXPECT_ERR_WARN(cbs, false, false);

        auto u = conf.getConfSubtree("u");
        EXPECT_ERR_WARN(cbs, false, false);

        EXPECT_EQ("data", u.getValue<std::string>());
        EXPECT_ERR_WARN(cbs, false, false);

        // test that read status of data is transferred in move construction
        {
            BaseLib::ConfigTreeNew u2(std::move(u));
            EXPECT_ERR_WARN(cbs, false, false);
        }
        EXPECT_ERR_WARN(cbs, false, false);

        // test that read status of children is transferred in move construction
        BaseLib::ConfigTreeNew conf2(std::move(conf));

        EXPECT_EQ("XX",   conf2.getConfParam<std::string>("n", "XX"));
        EXPECT_ERR_WARN(cbs, false, false);

        conf2.checkConfParam("t", "Test");
        EXPECT_ERR_WARN(cbs, false, false);

        BaseLib::checkAndInvalidate(conf2);
        EXPECT_ERR_WARN(cbs, false, false);
    } // ConfigTree destroyed here
    EXPECT_ERR_WARN(cbs, false, false);
}

// String literals are somewhat special for template classes
TEST(BaseLibConfigTree, MoveAssign)
{
    const char xml[] =
            "<s>test</s>"
            "<t>Test</t>"
            "<u>data</u>";
    auto const ptree = readXml(xml);

    Callbacks cbs;
    {
        BaseLib::ConfigTreeNew conf(ptree, cbs.get_error_cb(), cbs.get_warning_cb());

        EXPECT_EQ("test", conf.getConfParam<std::string>("s", "XX"));
        EXPECT_ERR_WARN(cbs, false, false);

        auto u = conf.getConfSubtree("u");
        EXPECT_ERR_WARN(cbs, false, false);

        EXPECT_EQ("data", u.getValue<std::string>());
        EXPECT_ERR_WARN(cbs, false, false);

        // test that read status of data is transferred in move assignment
        {
            BaseLib::ConfigTreeNew u2(ptree, cbs.get_error_cb(), cbs.get_warning_cb());
            u2 = std::move(u);
            // Expect warning because u2 has not been traversed
            // entirely before assignment.
            EXPECT_ERR_WARN(cbs, false, true);
        }
        EXPECT_ERR_WARN(cbs, false, false);

        // test that read status of children is transferred in move construction
        {
            BaseLib::ConfigTreeNew conf2(ptree, cbs.get_error_cb(), cbs.get_warning_cb());
            conf2 = std::move(conf);
            // Expect warning because conf2 has not been traversed
            // entirely before assignment.
            EXPECT_ERR_WARN(cbs, false, true);

            EXPECT_EQ("XX",   conf2.getConfParam<std::string>("n", "XX"));
            EXPECT_ERR_WARN(cbs, false, false);

            conf2.checkConfParam("t", "Test");
            EXPECT_ERR_WARN(cbs, false, false);
        }
        EXPECT_ERR_WARN(cbs, false, false);
    } // ConfigTree destroyed here
    EXPECT_ERR_WARN(cbs, false, false);
}
