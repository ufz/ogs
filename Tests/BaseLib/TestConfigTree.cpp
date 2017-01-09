/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>
#include <logog/include/logog.hpp>

#include <boost/property_tree/xml_parser.hpp>
#include <numeric>
#include <sstream>
#include <vector>

#include "BaseLib/ConfigTree.h"
#include "Tests/TestTools.h"

// make useful line numbers appear in the output of gtest
#define EXPECT_ERR_WARN(cbs, error, warning) do  { \
        if (error)   EXPECT_TRUE((cbs).get_error());   else EXPECT_FALSE((cbs).get_error()); \
        if (warning) EXPECT_TRUE((cbs).get_warning()); else EXPECT_FALSE((cbs).get_warning()); \
        (cbs).reset(); \
    } while(false)

// Exception thrown by the error callback of the class below
class Exc {};

// class that provides callback functions used with ConfigTree
class Callbacks
{
public:
    BaseLib::ConfigTree::Callback
    get_error_cb()
    {
        return [this](std::string const& filename, std::string const& path,
                      std::string const& message)
        {
            (void) path; (void) message;

            // check that filename is passed around properly, especially with
            // move construction/assignment
            EXPECT_EQ("FILENAME", filename);

            DBUG("error <%s> : %s", path.c_str(), message.c_str());
            _error = true;
            throw Exc(); // throw in order to stop normal execution
        };
    }

    BaseLib::ConfigTree::Callback
    get_warning_cb()
    {
        return [this](std::string const& filename, std::string const& path,
                      std::string const& message)
        {
            (void) path; (void) message;

            // check that filename is passed around properly, especially with
            // move construction/assignment
            EXPECT_EQ("FILENAME", filename);

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

BaseLib::ConfigTree
makeConfigTree(boost::property_tree::ptree const& ptree, Callbacks& cbs)
{
    return BaseLib::ConfigTree(ptree, "FILENAME",
                                  cbs.get_error_cb(), cbs.get_warning_cb());
}


TEST(BaseLibConfigTree, Empty)
{
    boost::property_tree::ptree ptree;
    Callbacks cbs;

    {
        auto const conf = makeConfigTree(ptree, cbs);
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
            "<z attr=\"0.5\" optattr=\"false\">32.0</z>"
            "<vector>0 1 2 3 4</vector>"
            "<vector_bad1>x 1 2a</vector_bad1>"
            "<vector_bad2>0 1 2a</vector_bad2>"
            ;
    auto const ptree = readXml(xml);

    Callbacks cbs;
    {
        auto const conf = makeConfigTree(ptree, cbs);

        EXPECT_EQ(5.6e-4, conf.getConfigParameter<double>("double")); // read certain types
        EXPECT_ERR_WARN(cbs, false, false);
        EXPECT_TRUE(conf.getConfigParameter<bool>("bool"));
        EXPECT_ERR_WARN(cbs, false, false);
        EXPECT_EQ(5, conf.getConfigParameter<int>("int"));
        EXPECT_ERR_WARN(cbs, false, false);

        EXPECT_EQ(8, conf.getConfigParameter<int>("intx", 8)); // reading with default value
        EXPECT_ERR_WARN(cbs, false, false);


        // Testing subtree
        {
            auto sub = conf.getConfigSubtree("sub");
            EXPECT_ERR_WARN(cbs, false, false);

            EXPECT_EQ(6.1f, sub.getConfigParameter<float>("float"));
            EXPECT_ERR_WARN(cbs, false, false);

            if (auto f2 = sub.getConfigParameterOptional<float>("float2")) { // read optional value
                EXPECT_EQ(0.1f, *f2);
            }
            EXPECT_ERR_WARN(cbs, false, false);

            auto f3 = sub.getConfigParameterOptional<float>("float3"); // optional value not existent
            ASSERT_FALSE(f3);
            EXPECT_ERR_WARN(cbs, false, false);


            // Testing the getConfigParameter...() (non-template) / getValue() combination

            auto bool1 = sub.getConfigSubtree("bool1");
            EXPECT_ERR_WARN(cbs, false, false);
            EXPECT_FALSE(bool1.getValue<bool>());
            EXPECT_ERR_WARN(cbs, false, false);
            EXPECT_ANY_THROW(bool1.getValue<bool>()); // getting data twice
            EXPECT_ERR_WARN(cbs, true, false);

            if (auto bool2 = sub.getConfigSubtreeOptional("bool2")) {
                EXPECT_ERR_WARN(cbs, false, false);
                EXPECT_FALSE(bool2->getValue<bool>());
            }
            EXPECT_ERR_WARN(cbs, false, false);

            if (auto bool3 = sub.getConfigSubtreeOptional("bool3")) {
                EXPECT_ERR_WARN(cbs, false, false);
                EXPECT_ANY_THROW(bool3->getValue<bool>());
                EXPECT_ERR_WARN(cbs, true, false); // error because of no data
            }
            EXPECT_ERR_WARN(cbs, false, false);

            EXPECT_FALSE(sub.getConfigSubtreeOptional("bool4")); // optional value not existent
            EXPECT_ERR_WARN(cbs, false, false);


            // Testing ignore

            sub.ignoreConfigParameter("ignored");
            EXPECT_ERR_WARN(cbs, false, false);
            sub.ignoreConfigParameterAll("ignored2");
            EXPECT_ERR_WARN(cbs, false, false);
            sub.ignoreConfigParameterAll("ignored4"); // I can ignore nonexistent stuff
            EXPECT_ERR_WARN(cbs, false, false);

            // I can not ignore stuff that I already read
            // this also makes sure that the subtree inherits the callbacks properly
            EXPECT_ANY_THROW(sub.ignoreConfigParameter("float"));
            EXPECT_ERR_WARN(cbs, true, false);
        }
        for (int i : {0, 1, 2}) {
            (void) i;
            EXPECT_EQ("Y", conf.peekConfigParameter<std::string>("x"));
            EXPECT_ERR_WARN(cbs, false, false);
        }
        conf.checkConfigParameter<std::string>("x", "Y");
        EXPECT_ERR_WARN(cbs, false, false);


        // Testing attributes
        {
            auto z = conf.getConfigSubtree("z");
            EXPECT_ERR_WARN(cbs, false, false);
            EXPECT_EQ(0.5, z.getConfigAttribute<double>("attr"));
            EXPECT_ERR_WARN(cbs, false, false);
            EXPECT_ANY_THROW(z.getConfigAttribute<double>("attr")); // getting attribute twice
            EXPECT_ERR_WARN(cbs, true, false);
            EXPECT_ANY_THROW(z.getConfigAttribute<double>("not_an_attr")); // nonexistent attribute
            EXPECT_ERR_WARN(cbs, true, false);
            EXPECT_EQ(32.0, z.getValue<double>());
            EXPECT_ERR_WARN(cbs, false, false);
            auto const opt = z.getConfigAttributeOptional<bool>("optattr");
            EXPECT_TRUE(!!opt); EXPECT_FALSE(*opt);
            EXPECT_ERR_WARN(cbs, false, false);
            EXPECT_ANY_THROW(z.getConfigAttributeOptional<bool>("optattr")); // getting attribute twice
            EXPECT_ERR_WARN(cbs, true, false);
            EXPECT_FALSE(z.getConfigAttributeOptional<bool>("also_not_an_attr")); // nonexisting attribute
            EXPECT_ERR_WARN(cbs, false, false);
        }

        // Testing vector
        {
            auto v = conf.getConfigParameter<std::vector<int>>("vector");
            EXPECT_ERR_WARN(cbs, false, false);
            EXPECT_EQ(5u, v.size());
            std::vector<int> expected_vector(5);
            std::iota(expected_vector.begin(), expected_vector.end(), 0);
            EXPECT_TRUE(std::equal(expected_vector.begin(),
                                   expected_vector.end(), v.begin()));
            EXPECT_ANY_THROW(
                conf.getConfigParameter<std::vector<int>>("vector_bad1"));
            EXPECT_ERR_WARN(cbs, true, false);
            EXPECT_ANY_THROW(
                conf.getConfigParameter<std::vector<int>>("vector_bad2"));
            EXPECT_ERR_WARN(cbs, true, false);
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
        auto const conf = makeConfigTree(ptree, cbs);

        EXPECT_EQ(5.6, conf.getConfigParameter<double>("double"));
        EXPECT_ERR_WARN(cbs, false, false);

        conf.getConfigSubtree("tag");
        EXPECT_ERR_WARN(cbs, false, true); // data of <tag> has not been read

        EXPECT_EQ(1, conf.getConfigParameter<int>("pt"));
        EXPECT_ERR_WARN(cbs, false, true); // attribute "x" has not been read

        {
            auto pt2 = conf.getConfigSubtree("pt2");
            EXPECT_EQ(0.5, pt2.getConfigAttribute<double>("x"));
            EXPECT_ERR_WARN(cbs, false, false);
            EXPECT_EQ(1.0, pt2.getConfigAttribute<double>("y"));
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
        auto const conf = makeConfigTree(ptree, cbs);

        {
            // check that std::distance can be computed twice in a row
            auto list = conf.getConfigSubtreeList("val");
            EXPECT_ERR_WARN(cbs, false, false);
            EXPECT_EQ(3, std::distance(list.begin(), list.end()));
            EXPECT_ERR_WARN(cbs, false, false);
            EXPECT_EQ(3, std::distance(list.begin(), list.end()));
            EXPECT_ERR_WARN(cbs, false, false);
        }

        {
            // check that std::distance can be computed twice in a row
            auto list = conf.getConfigParameterList<int>("int");
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
        auto const conf = makeConfigTree(ptree, cbs);

        for (auto p : conf.getConfigSubtreeList("nonexistent_list"))
        {
            (void) p;
            FAIL() << "Expected empty list";
        }
        EXPECT_ERR_WARN(cbs, false, false);

        int i = 0;
        for (auto ct : conf.getConfigSubtreeList("val"))
        {
            EXPECT_EQ(i, ct.getConfigParameter<int>("int"));
            EXPECT_ERR_WARN(cbs, false, false);
            ++i;
        }
    } // ConfigTree destroyed here
    EXPECT_ERR_WARN(cbs, false, false);
}

TEST(BaseLibConfigTree, GetParamList)
{
    const char xml[] =
            "<int>0</int>"
            "<int>1</int>"
            "<int>2</int>"
            "<int2 a=\"b\">3</int2>"
            "<int3>4<error/></int3>";
    auto const ptree = readXml(xml);

    Callbacks cbs;
    {
        auto const conf = makeConfigTree(ptree, cbs);

        for (auto p : conf.getConfigParameterList("nonexistent_list"))
        {
            (void) p;
            FAIL() << "Expected empty list";
        }
        EXPECT_ERR_WARN(cbs, false, false);

        int i = 0;
        for (auto p : conf.getConfigParameterList("int"))
        {
            EXPECT_EQ(i, p.getValue<int>());
            EXPECT_ERR_WARN(cbs, false, false);
            ++i;
        }

        for (auto p : conf.getConfigParameterList("int2"))
        {
            EXPECT_EQ(i, p.getValue<int>());
            EXPECT_ERR_WARN(cbs, false, false);
            ++i;
        }
        EXPECT_ERR_WARN(cbs, false, true); // attribute "a" not read

        {
            // get list of parameters, i.e., subtrees without children
            auto range = conf.getConfigParameterList("int3");
            EXPECT_ERR_WARN(cbs, false, false);

            EXPECT_ANY_THROW(*range.begin());
            // error because of child tag <error/>
            EXPECT_ERR_WARN(cbs, true, false);
        } // range destroyed here
        EXPECT_ERR_WARN(cbs, false, false);

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
        auto const conf = makeConfigTree(ptree, cbs);

        for (auto p : conf.getConfigParameterList<int>("nonexistent_list"))
        {
            (void) p;
            FAIL() << "Expected empty list";
        }
        EXPECT_ERR_WARN(cbs, false, false);

        int n = 0;
        for (auto i : conf.getConfigParameterList<int>("int"))
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
        auto const conf = makeConfigTree(ptree, cbs);

        EXPECT_ANY_THROW(conf.getConfigParameter<int>("int"));
        EXPECT_ERR_WARN(cbs, true, false);
        EXPECT_ANY_THROW(conf.ignoreConfigParameter("int")); // after failure I also cannot ignore something
        EXPECT_ERR_WARN(cbs, true, false);

        EXPECT_ANY_THROW(conf.getConfigParameter<double>("double"));
        EXPECT_ERR_WARN(cbs, true, false);

        // peek value existent but not convertible
        EXPECT_ANY_THROW(conf.peekConfigParameter<double>("non_double"));
        EXPECT_ERR_WARN(cbs, true, false);

        // optional value existent but not convertible
        EXPECT_ANY_THROW(
            auto d = conf.getConfigParameterOptional<double>("non_double");
            ASSERT_FALSE(d);
        );
        EXPECT_ERR_WARN(cbs, true, false);

        // assert that I can only ignore something once
        conf.ignoreConfigParameter("ign");
        EXPECT_ERR_WARN(cbs, false, false);
        EXPECT_ANY_THROW(conf.ignoreConfigParameter("ign"));
        EXPECT_ERR_WARN(cbs, true, false);
        conf.ignoreConfigParameterAll("ign2");
        EXPECT_ERR_WARN(cbs, false, false);
        EXPECT_ANY_THROW(conf.ignoreConfigParameterAll("ign2"));
        EXPECT_ERR_WARN(cbs, true, false);

        // assert that I cannot read a parameter twice
        conf.getConfigParameter<bool>("bool");
        EXPECT_ERR_WARN(cbs, false, false);
        EXPECT_ANY_THROW(conf.getConfigParameter<bool>("bool"));
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
        auto const conf = makeConfigTree(ptree, cbs);

        for (auto tag : { "<", "Z", ".", "$", "0", "", "/", "_", "a__" })
        {
            EXPECT_ANY_THROW(conf.getConfigParameter<int>(tag));
            EXPECT_ERR_WARN(cbs, true, false);
            EXPECT_ANY_THROW(conf.getConfigParameter<int>(tag, 500));
            EXPECT_ERR_WARN(cbs, true, false);
            EXPECT_ANY_THROW(conf.getConfigParameterOptional<int>(tag));
            EXPECT_ERR_WARN(cbs, true, false);
            EXPECT_ANY_THROW(conf.getConfigParameterList<int>(tag));
            EXPECT_ERR_WARN(cbs, true, false);

            EXPECT_ANY_THROW(conf.peekConfigParameter<int>(tag));
            EXPECT_ERR_WARN(cbs, true, false);
            EXPECT_ANY_THROW(conf.checkConfigParameter<int>(tag, 500));
            EXPECT_ERR_WARN(cbs, true, false);

            EXPECT_ANY_THROW(conf.getConfigSubtree(tag));
            EXPECT_ERR_WARN(cbs, true, false);
            EXPECT_ANY_THROW(conf.getConfigSubtreeOptional(tag));
            EXPECT_ERR_WARN(cbs, true, false);
            EXPECT_ANY_THROW(conf.getConfigSubtreeList(tag));
            EXPECT_ERR_WARN(cbs, true, false);

            EXPECT_ANY_THROW(conf.getConfigAttribute<int>(tag));
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
        auto const conf = makeConfigTree(ptree, cbs);

        EXPECT_EQ("test", conf.getConfigParameter<std::string>("s", "XX"));
        EXPECT_ERR_WARN(cbs, false, false);

        // <n> not present in the XML, so return the default value
        EXPECT_EQ("XX",   conf.getConfigParameter<std::string>("n", "XX"));
        EXPECT_ERR_WARN(cbs, false, false);

        conf.checkConfigParameter("t", "Test");
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
        auto conf = makeConfigTree(ptree, cbs);

        EXPECT_EQ("test", conf.getConfigParameter<std::string>("s", "XX"));
        EXPECT_ERR_WARN(cbs, false, false);

        auto u = conf.getConfigSubtree("u");
        EXPECT_ERR_WARN(cbs, false, false);

        EXPECT_EQ("data", u.getValue<std::string>());
        EXPECT_ERR_WARN(cbs, false, false);

        // test that read status of data is transferred in move construction
        {
            BaseLib::ConfigTree const u2(std::move(u));
            EXPECT_ERR_WARN(cbs, false, false);
        }
        EXPECT_ERR_WARN(cbs, false, false);

        // test that read status of children is transferred in move construction
        BaseLib::ConfigTree conf2(std::move(conf));

        EXPECT_EQ("XX",   conf2.getConfigParameter<std::string>("n", "XX"));
        EXPECT_ERR_WARN(cbs, false, false);

        conf2.checkConfigParameter("t", "Test");
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
        auto conf = makeConfigTree(ptree, cbs);

        EXPECT_EQ("test", conf.getConfigParameter<std::string>("s", "XX"));
        EXPECT_ERR_WARN(cbs, false, false);

        auto u = conf.getConfigSubtree("u");
        EXPECT_ERR_WARN(cbs, false, false);

        EXPECT_EQ("data", u.getValue<std::string>());
        EXPECT_ERR_WARN(cbs, false, false);

        // test that read status of data is transferred in move assignment
        {
            auto u2 = makeConfigTree(ptree, cbs);
            u2 = std::move(u);
            // Expect warning because u2 has not been traversed
            // entirely before assignment.
            EXPECT_ERR_WARN(cbs, false, true);
        }
        EXPECT_ERR_WARN(cbs, false, false);

        // test that read status of children is transferred in move construction
        {
            auto conf2 = makeConfigTree(ptree, cbs);
            conf2 = std::move(conf);
            // Expect warning because conf2 has not been traversed
            // entirely before assignment.
            EXPECT_ERR_WARN(cbs, false, true);

            EXPECT_EQ("XX",   conf2.getConfigParameter<std::string>("n", "XX"));
            EXPECT_ERR_WARN(cbs, false, false);

            conf2.checkConfigParameter("t", "Test");
            EXPECT_ERR_WARN(cbs, false, false);
        }
        EXPECT_ERR_WARN(cbs, false, false);
    } // ConfigTree destroyed here
    EXPECT_ERR_WARN(cbs, false, false);
}
