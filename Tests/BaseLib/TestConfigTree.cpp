// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <boost/property_tree/xml_parser.hpp>
#include <numeric>
#include <sstream>
#include <vector>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Logging.h"
#include "Tests/TestTools.h"

// make useful line numbers appear in the output of gtest
#define EXPECT_ERR_WARN(cbs, error, warning)   \
    do                                         \
    {                                          \
        if (error)                             \
            EXPECT_TRUE((cbs).get_error());    \
        else                                   \
            EXPECT_FALSE((cbs).get_error());   \
        if (warning)                           \
            EXPECT_TRUE((cbs).get_warning());  \
        else                                   \
            EXPECT_FALSE((cbs).get_warning()); \
        (cbs).reset();                         \
    } while (false)

// Exception thrown by the error callback of the class below
class Exc
{
};

// class that provides callback functions used with ConfigTree
class Callbacks
{
public:
    BaseLib::ConfigTree::Callback get_error_cb()
    {
        return [this](std::string const& filename, std::string const& path,
                      std::string const& message)
        {
            (void)path;
            (void)message;

            // check that filename is passed around properly, especially with
            // move construction/assignment
            EXPECT_EQ(std::filesystem::current_path() / "FILENAME", filename);

            DBUG("error <{:s}> : {:s}", path, message);
            _error = true;
            throw Exc();  // throw in order to stop normal execution
        };
    }

    BaseLib::ConfigTree::Callback get_warning_cb()
    {
        return [this](std::string const& filename, std::string const& path,
                      std::string const& message)
        {
            (void)path;
            (void)message;

            // check that filename is passed around properly, especially with
            // move construction/assignment
            EXPECT_EQ(std::filesystem::current_path() / "FILENAME", filename);

            DBUG("warning <{:s}> : {:s}", path, message);
            _warning = true;
        };
    }

    bool get_error() const { return _error; }
    bool get_warning() const { return _warning; }
    void reset()
    {
        _error = false;
        _warning = false;
    }

private:
    bool _error = false;
    bool _warning = false;
};

BaseLib::ConfigTree makeConfigTree(boost::property_tree::ptree&& ptree,
                                   Callbacks& cbs)
{
    return BaseLib::ConfigTree(std::move(ptree), "FILENAME", cbs.get_error_cb(),
                               cbs.get_warning_cb());
}

TEST(BaseLibConfigTree, Empty)
{
    boost::property_tree::ptree ptree;
    Callbacks cbs;

    {
        auto const conf = makeConfigTree(std::move(ptree), cbs);
        (void)conf;
    }  // ConfigTree destroyed here

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
        "<vector_bad2>0 1 2a</vector_bad2>";
    auto ptree = Tests::readXml(xml);

    Callbacks cbs;
    {
        auto const conf = makeConfigTree(std::move(ptree), cbs);

        EXPECT_EQ(5.6e-4, conf.getConfigParameter<double>(
                              "double"));  // read certain types
        EXPECT_ERR_WARN(cbs, false, false);
        EXPECT_TRUE(conf.getConfigParameter<bool>("bool"));
        EXPECT_ERR_WARN(cbs, false, false);
        EXPECT_EQ(5, conf.getConfigParameter<int>("int"));
        EXPECT_ERR_WARN(cbs, false, false);

        EXPECT_EQ(8, conf.getConfigParameter<int>(
                         "intx", 8));  // reading with default value
        EXPECT_ERR_WARN(cbs, false, false);

        // Testing subtree
        {
            auto sub = conf.getConfigSubtree("sub");
            EXPECT_ERR_WARN(cbs, false, false);

            EXPECT_EQ(6.1f, sub.getConfigParameter<float>("float"));
            EXPECT_ERR_WARN(cbs, false, false);

            if (auto f2 = sub.getConfigParameterOptional<float>("float2"))
            {  // read optional value
                EXPECT_EQ(0.1f, *f2);
            }
            EXPECT_ERR_WARN(cbs, false, false);

            auto f3 = sub.getConfigParameterOptional<float>(
                "float3");  // optional value not existent
            ASSERT_FALSE(f3);
            EXPECT_ERR_WARN(cbs, false, false);

            // Testing the getConfigParameter...() (non-template) / getValue()
            // combination

            auto bool1 = sub.getConfigSubtree("bool1");
            EXPECT_ERR_WARN(cbs, false, false);
            EXPECT_FALSE(bool1.getValue<bool>());
            EXPECT_ERR_WARN(cbs, false, false);
            EXPECT_ANY_THROW(bool1.getValue<bool>());  // getting data twice
            EXPECT_ERR_WARN(cbs, true, false);

            if (auto bool2 = sub.getConfigSubtreeOptional("bool2"))
            {
                EXPECT_ERR_WARN(cbs, false, false);
                EXPECT_FALSE(bool2->getValue<bool>());
            }
            EXPECT_ERR_WARN(cbs, false, false);

            if (auto bool3 = sub.getConfigSubtreeOptional("bool3"))
            {
                EXPECT_ERR_WARN(cbs, false, false);
                EXPECT_ANY_THROW(bool3->getValue<bool>());
                EXPECT_ERR_WARN(cbs, true, false);  // error because of no data
            }
            EXPECT_ERR_WARN(cbs, false, false);

            EXPECT_FALSE(sub.getConfigSubtreeOptional(
                "bool4"));  // optional value not existent
            EXPECT_ERR_WARN(cbs, false, false);

            // Testing ignore

            sub.ignoreConfigParameter("ignored");
            EXPECT_ERR_WARN(cbs, false, false);
            sub.ignoreConfigParameterAll("ignored2");
            EXPECT_ERR_WARN(cbs, false, false);
            sub.ignoreConfigParameterAll(
                "ignored4");  // I can ignore nonexistent stuff
            EXPECT_ERR_WARN(cbs, false, false);

            // I can not ignore stuff that I already read
            // this also makes sure that the subtree inherits the callbacks
            // properly
            EXPECT_ANY_THROW(sub.ignoreConfigParameter("float"));
            EXPECT_ERR_WARN(cbs, true, false);
        }
        for (int i : {0, 1, 2})
        {
            (void)i;
            EXPECT_EQ("Y", conf.peekConfigParameter<std::string>("x"));
            EXPECT_ERR_WARN(cbs, false, false);
        }
        conf.checkConfigParameter("x", "Y");
        EXPECT_ERR_WARN(cbs, false, false);

        // Testing attributes
        {
            auto z = conf.getConfigSubtree("z");
            EXPECT_ERR_WARN(cbs, false, false);
            EXPECT_EQ(0.5, z.getConfigAttribute<double>("attr"));
            EXPECT_ERR_WARN(cbs, false, false);
            EXPECT_ANY_THROW(z.getConfigAttribute<double>(
                "attr"));  // getting attribute twice
            EXPECT_ERR_WARN(cbs, true, false);
            EXPECT_ANY_THROW(z.getConfigAttribute<double>(
                "not_an_attr"));  // nonexistent attribute
            EXPECT_ERR_WARN(cbs, true, false);
            EXPECT_EQ(32.0, z.getValue<double>());
            EXPECT_ERR_WARN(cbs, false, false);
            auto const opt = z.getConfigAttributeOptional<bool>("optattr");
            EXPECT_TRUE(!!opt);
            EXPECT_FALSE(*opt);
            EXPECT_ERR_WARN(cbs, false, false);
            EXPECT_ANY_THROW(z.getConfigAttributeOptional<bool>(
                "optattr"));  // getting attribute twice
            EXPECT_ERR_WARN(cbs, true, false);
            EXPECT_FALSE(z.getConfigAttributeOptional<bool>(
                "also_not_an_attr"));  // nonexisting attribute
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
    }  // ConfigTree destroyed here
    EXPECT_ERR_WARN(cbs, false, false);
}

TEST(BaseLibConfigTree, IncompleteParse)
{
    const char xml[] =
        "<double>5.6</double>"
        "<not_read>true</not_read>"
        "<tag>this data won't be read</tag>"
        "<pt x=\"0.5\">1</pt>"
        "<pt2 x=\"0.5\" y=\"1.0\" z=\"2.0\" />";
    auto ptree = Tests::readXml(xml);

    Callbacks cbs;
    {
        auto const conf = makeConfigTree(std::move(ptree), cbs);

        EXPECT_EQ(5.6, conf.getConfigParameter<double>("double"));
        EXPECT_ERR_WARN(cbs, false, false);

        conf.getConfigSubtree("tag");
        EXPECT_ERR_WARN(cbs, false, true);  // data of <tag> has not been read

        EXPECT_EQ(1, conf.getConfigParameter<int>("pt"));
        EXPECT_ERR_WARN(cbs, false, true);  // attribute "x" has not been read

        {
            auto pt2 = conf.getConfigSubtree("pt2");
            EXPECT_EQ(0.5, pt2.getConfigAttribute<double>("x"));
            EXPECT_ERR_WARN(cbs, false, false);
            EXPECT_EQ(1.0, pt2.getConfigAttribute<double>("y"));
            EXPECT_ERR_WARN(cbs, false, false);

            BaseLib::checkAndInvalidate(pt2);
            EXPECT_ERR_WARN(cbs, false, true);  // attribute "z" not read
        }
        EXPECT_ERR_WARN(cbs, false, false);

    }  // ConfigTree destroyed here
    EXPECT_ERR_WARN(cbs, false,
                    true);  // expect warning because I didn't read everything
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
    auto ptree = Tests::readXml(xml);

    Callbacks cbs;
    {
        auto const conf = makeConfigTree(std::move(ptree), cbs);

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

    }  // ConfigTree destroyed here

    // there will be warnings because I don't process the list entries
    EXPECT_ERR_WARN(cbs, false, true);
}

TEST(BaseLibConfigTree, GetSubtreeList)
{
    const char xml[] =
        "<val><int>0</int></val>"
        "<val><int>1</int></val>"
        "<val><int>2</int></val>";
    auto ptree = Tests::readXml(xml);

    Callbacks cbs;
    {
        auto const conf = makeConfigTree(std::move(ptree), cbs);

        auto const expected_empty_list =
            conf.getConfigParameterList("nonexistent_list");
        ASSERT_TRUE(expected_empty_list.empty()) << "Expected empty list";

        EXPECT_ERR_WARN(cbs, false, false);

        int i = 0;
        for (auto ct : conf.getConfigSubtreeList("val"))
        {
            EXPECT_EQ(i, ct.getConfigParameter<int>("int"));
            EXPECT_ERR_WARN(cbs, false, false);
            ++i;
        }
    }  // ConfigTree destroyed here
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
    auto ptree = Tests::readXml(xml);

    Callbacks cbs;
    {
        auto const conf = makeConfigTree(std::move(ptree), cbs);

        auto const expected_empty_list =
            conf.getConfigParameterList("nonexistent_list");
        ASSERT_TRUE(expected_empty_list.empty()) << "Expected empty list";

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
        EXPECT_ERR_WARN(cbs, false, true);  // attribute "a" not read

        {
            // get list of parameters, i.e., subtrees without children
            auto range = conf.getConfigParameterList("int3");
            EXPECT_ERR_WARN(cbs, false, false);

            EXPECT_ANY_THROW(*range.begin());
            // error because of child tag <error/>
            EXPECT_ERR_WARN(cbs, true, false);
        }  // range destroyed here
        EXPECT_ERR_WARN(cbs, false, false);

    }  // ConfigTree destroyed here
    EXPECT_ERR_WARN(cbs, false, false);
}

TEST(BaseLibConfigTree, GetValueList)
{
    const char xml[] =
        "<int>0</int>"
        "<int>1</int>"
        "<int>2</int>";
    auto ptree = Tests::readXml(xml);

    Callbacks cbs;
    {
        auto const conf = makeConfigTree(std::move(ptree), cbs);

        auto const expected_empty_list =
            conf.getConfigParameterList<int>("nonexistent_list");
        ASSERT_TRUE(expected_empty_list.empty()) << "Expected empty list";

        EXPECT_ERR_WARN(cbs, false, false);

        int n = 0;
        for (auto i : conf.getConfigParameterList<int>("int"))
        {
            EXPECT_EQ(n, i);
            EXPECT_ERR_WARN(cbs, false, false);
            ++n;
        }
    }  // ConfigTree destroyed here
    EXPECT_ERR_WARN(cbs, false, false);
}

TEST(BaseLibConfigTree, NoConversion)
{
    const char xml[] =
        "<int>5.6</int>"                 // not convertible to int
        "<double>5.6tz</double>"         // not convertible to double
        "<non_double>0.1x</non_double>"  // not either convertible to double
        "<bool>true</bool>"
        "<ign/>"
        "<ign2/><ign2/><ign2/>";
    auto ptree = Tests::readXml(xml);

    Callbacks cbs;
    {
        auto const conf = makeConfigTree(std::move(ptree), cbs);

        EXPECT_ANY_THROW(conf.getConfigParameter<int>("int"));
        EXPECT_ERR_WARN(cbs, true, false);
        EXPECT_ANY_THROW(conf.ignoreConfigParameter(
            "int"));  // after failure I also cannot ignore something
        EXPECT_ERR_WARN(cbs, true, false);

        EXPECT_ANY_THROW(conf.getConfigParameter<double>("double"));
        EXPECT_ERR_WARN(cbs, true, false);

        // peek value existent but not convertible
        EXPECT_ANY_THROW(conf.peekConfigParameter<double>("non_double"));
        EXPECT_ERR_WARN(cbs, true, false);

        // optional value existent but not convertible
        EXPECT_ANY_THROW(
            auto d = conf.getConfigParameterOptional<double>("non_double");
            ASSERT_FALSE(d););
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

    }  // ConfigTree destroyed here

    // There will bewarnings because I don't succeed in reading every setting,
    // and furthermore I read some setting too often.
    EXPECT_ERR_WARN(cbs, false, false);
}

TEST(BaseLibConfigTree, BadKeynames)
{
    const char xml[] = "";
    auto ptree = Tests::readXml(xml);

    Callbacks cbs;
    {
        auto const conf = makeConfigTree(std::move(ptree), cbs);

        for (auto tag : {"<", "Z", ".", "$", "0", "", "/", "_", "a__"})
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
            EXPECT_ANY_THROW(conf.checkConfigParameter(tag, "500"));
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

    }  // ConfigTree destroyed here

    EXPECT_ERR_WARN(cbs, false, false);
}

// String literals are somewhat special for template classes
TEST(BaseLibConfigTree, StringLiterals)
{
    const char xml[] =
        "<s>test</s>"
        "<t>Test</t>";
    auto ptree = Tests::readXml(xml);

    Callbacks cbs;
    {
        auto const conf = makeConfigTree(std::move(ptree), cbs);

        EXPECT_EQ("test", conf.getConfigParameter<std::string>("s", "XX"));
        EXPECT_ERR_WARN(cbs, false, false);

        // <n> not present in the XML, so return the default value
        EXPECT_EQ("XX", conf.getConfigParameter<std::string>("n", "XX"));
        EXPECT_ERR_WARN(cbs, false, false);

        conf.checkConfigParameter("t", "Test");
        EXPECT_ERR_WARN(cbs, false, false);
    }  // ConfigTree destroyed here
    EXPECT_ERR_WARN(cbs, false, false);
}

// String literals are somewhat special for template classes
TEST(BaseLibConfigTree, MoveConstruct)
{
    const char xml[] =
        "<s>test</s>"
        "<t>Test</t>"
        "<u>data</u>";
    auto ptree = Tests::readXml(xml);

    Callbacks cbs;
    {
        auto conf = makeConfigTree(std::move(ptree), cbs);

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

        EXPECT_EQ("XX", conf2.getConfigParameter<std::string>("n", "XX"));
        EXPECT_ERR_WARN(cbs, false, false);

        conf2.checkConfigParameter("t", "Test");
        EXPECT_ERR_WARN(cbs, false, false);

        BaseLib::checkAndInvalidate(conf2);
        EXPECT_ERR_WARN(cbs, false, false);
    }  // ConfigTree destroyed here
    EXPECT_ERR_WARN(cbs, false, false);
}

// String literals are somewhat special for template classes
TEST(BaseLibConfigTree, MoveAssign)
{
    const char xml[] =
        "<s>test</s>"
        "<t>Test</t>"
        "<u>data</u>";
    auto const ptree = Tests::readXml(xml);

    Callbacks cbs;
    {
        auto conf = makeConfigTree(boost::property_tree::ptree(ptree), cbs);

        EXPECT_EQ("test", conf.getConfigParameter<std::string>("s", "XX"));
        EXPECT_ERR_WARN(cbs, false, false);

        auto u = conf.getConfigSubtree("u");
        EXPECT_ERR_WARN(cbs, false, false);

        EXPECT_EQ("data", u.getValue<std::string>());
        EXPECT_ERR_WARN(cbs, false, false);

        // test that read status of data is transferred in move assignment
        {
            auto u2 = makeConfigTree(boost::property_tree::ptree(ptree), cbs);
            u2 = std::move(u);
            // Expect warning because u2 has not been traversed
            // entirely before assignment.
            EXPECT_ERR_WARN(cbs, false, true);
        }
        EXPECT_ERR_WARN(cbs, false, false);

        // test that read status of children is transferred in move construction
        {
            auto conf2 =
                makeConfigTree(boost::property_tree::ptree(ptree), cbs);
            conf2 = std::move(conf);
            // Expect warning because conf2 has not been traversed
            // entirely before assignment.
            EXPECT_ERR_WARN(cbs, false, true);

            EXPECT_EQ("XX", conf2.getConfigParameter<std::string>("n", "XX"));
            EXPECT_ERR_WARN(cbs, false, false);

            conf2.checkConfigParameter("t", "Test");
            EXPECT_ERR_WARN(cbs, false, false);
        }
        EXPECT_ERR_WARN(cbs, false, false);
    }  // ConfigTree destroyed here
    EXPECT_ERR_WARN(cbs, false, false);
}

TEST(BaseLibConfigTree, GetAllChildrenEmptyTreeReturnsNoChildren)
{
    boost::property_tree::ptree ptree;
    Callbacks cbs;
    {
        auto const conf = makeConfigTree(std::move(ptree), cbs);
        auto children = conf.getAllChildren();
        EXPECT_TRUE(children.empty());
        EXPECT_ERR_WARN(cbs, false, false);
    }
    EXPECT_ERR_WARN(cbs, false, false);
}

// Every non-attribute child is returned, in document order; <xmlattr> storage
// is never reported as a child; and reading each child's full content
// (immediate data plus attributes) consumes it, so nothing warns.
TEST(BaseLibConfigTree, GetAllChildrenReturnsEveryNonAttributeChildInOrder)
{
    const char xml[] =
        "<val a=\"x\" b=\"y\">z</val>"
        "<val2>w</val2>"
        "<val3 p=\"q\"/>";
    auto ptree = Tests::readXml(xml);
    Callbacks cbs;
    {
        auto const conf = makeConfigTree(std::move(ptree), cbs);
        auto children = conf.getAllChildren();

        EXPECT_EQ(3, children.size());
        // Tag order must be preserved.
        EXPECT_EQ("val", children[0].first);
        EXPECT_EQ("val2", children[1].first);
        EXPECT_EQ("val3", children[2].first);

        // <xmlattr> must not be among the children.
        for (auto const& [tag, _] : children)
        {
            EXPECT_NE("<xmlattr>", tag);
        }

        // Consume every child's content: 'val' has data and two attributes,
        // 'val2' has data, 'val3' has only one attribute (and is still
        // returned as a child).
        EXPECT_EQ("z", children[0].second->getValue<std::string>());
        EXPECT_EQ("x",
                  children[0].second->getConfigAttribute<std::string>("a"));
        EXPECT_EQ("y",
                  children[0].second->getConfigAttribute<std::string>("b"));
        // Dereference via operator* (equivalent to operator->).
        EXPECT_EQ("w", (*children[1].second).getValue<std::string>());
        EXPECT_EQ("q",
                  children[2].second->getConfigAttribute<std::string>("p"));
        EXPECT_ERR_WARN(cbs, false, false);
    }
    EXPECT_ERR_WARN(cbs, false, false);
}

// Regression: peeking children without consuming them must warn, exactly like
// an unread getConfigSubtree() result. getAllChildren() hands out each child as
// its own subtree; a child left unread reports its own unread content on
// destruction.
TEST(BaseLibConfigTree, GetAllChildrenUnreadChildWarns)
{
    const char xml[] =
        "<a>1</a>"
        "<b>2</b>";
    auto ptree = Tests::readXml(xml);
    Callbacks cbs;
    {
        auto const conf = makeConfigTree(std::move(ptree), cbs);
        auto children = conf.getAllChildren();
        // Only inspect the tags; read none of the children's values.
        EXPECT_EQ(2, children.size());
        EXPECT_EQ("a", children[0].first);
        EXPECT_EQ("b", children[1].first);
        EXPECT_ERR_WARN(cbs, false, false);
    }
    // Both children carry unread immediate data -> warning.
    EXPECT_ERR_WARN(cbs, false, true);
}

// A child can be descended into: reading a nested subtree through the child
// (child->getConfigSubtree(...)) consumes it.
TEST(BaseLibConfigTree, GetAllChildrenChildCanBeDescendedInto)
{
    const char xml[] =
        "<outer1><inner>alpha</inner></outer1>"
        "<outer2><inner>beta</inner></outer2>";
    auto ptree = Tests::readXml(xml);
    Callbacks cbs;
    {
        auto const conf = makeConfigTree(std::move(ptree), cbs);
        auto children = conf.getAllChildren();

        EXPECT_EQ(2, children.size());
        auto inner1 = children[0].second->getConfigSubtree("inner");
        EXPECT_EQ("alpha", inner1.getValue<std::string>());
        auto inner2 = children[1].second->getConfigSubtree("inner");
        EXPECT_EQ("beta", inner2.getValue<std::string>());
        EXPECT_ERR_WARN(cbs, false, false);
    }
    EXPECT_ERR_WARN(cbs, false, false);
}

// getAllChildren() works on a subtree (not only the root), and duplicate
// sibling tags are each returned as a separate entry. The subtree's own
// <xmlattr> storage is skipped; consuming everything warns nothing.
TEST(BaseLibConfigTree, GetAllChildrenOnSubtreeWithDuplicateSiblingTags)
{
    const char xml[] =
        "<root a=\"1\">"
        "<dup>x</dup>"
        "<dup>y</dup>"
        "</root>";
    auto ptree = Tests::readXml(xml);
    Callbacks cbs;
    {
        auto const conf = makeConfigTree(std::move(ptree), cbs);
        auto const sub = conf.getConfigSubtree("root");
        // Consume the subtree's own attribute so it does not warn.
        EXPECT_EQ(1, sub.getConfigAttribute<int>("a"));

        auto children = sub.getAllChildren();
        // Only the two <dup> elements are returned; <xmlattr> is skipped.
        EXPECT_EQ(2, children.size());
        EXPECT_EQ("dup", children[0].first);
        EXPECT_EQ("x", children[0].second->getValue<std::string>());
        EXPECT_EQ("dup", children[1].first);
        EXPECT_EQ("y", children[1].second->getValue<std::string>());

        for (auto const& [tag, _] : children)
        {
            EXPECT_NE("<xmlattr>", tag);
        }
        EXPECT_ERR_WARN(cbs, false, false);
    }
    EXPECT_ERR_WARN(cbs, false, false);
}

// Partial consume: reading some children but not all leaves the unread ones to
// warn, while the consumed ones stay silent.
TEST(BaseLibConfigTree, GetAllChildrenPartialConsumeWarnsForUnread)
{
    const char xml[] =
        "<a>1</a>"
        "<b>2</b>"
        "<c>3</c>"
        "<d>4</d>";
    auto ptree = Tests::readXml(xml);
    Callbacks cbs;
    {
        auto const conf = makeConfigTree(std::move(ptree), cbs);
        auto children = conf.getAllChildren();
        // Consume 'a' and 'c'; 'b' and 'd' remain unconsumed.
        EXPECT_EQ(4, children.size());
        EXPECT_EQ("1", children[0].second->getValue<std::string>());
        EXPECT_EQ("3", children[2].second->getValue<std::string>());
        EXPECT_ERR_WARN(cbs, false, false);
    }
    // 'b' and 'd' are unconsumed children -> warning.
    EXPECT_ERR_WARN(cbs, false, true);
}

// Duplicate sibling tags with partial consume: each occurrence is tracked
// separately, so consuming one of three same-tag children leaves the other two
// unread -> warning.
TEST(BaseLibConfigTree, GetAllChildrenDuplicateTagsPartialConsumeWarns)
{
    const char xml[] =
        "<a>1</a>"
        "<a>2</a>"
        "<a>3</a>";
    auto ptree = Tests::readXml(xml);
    Callbacks cbs;
    {
        auto const conf = makeConfigTree(std::move(ptree), cbs);
        auto children = conf.getAllChildren();
        EXPECT_EQ(3, children.size());
        // Consume only the first <a>; the other two remain unread.
        EXPECT_EQ("1", children[0].second->getValue<std::string>());
        EXPECT_ERR_WARN(cbs, false, false);
    }
    // Two unread <a> siblings -> warning.
    EXPECT_ERR_WARN(cbs, false, true);
}

// A consumed child must consume ALL of its attributes: reading 'a' but not 'b'
// leaves 'b' unread -> warning.
TEST(BaseLibConfigTree, GetAllChildrenUnreadAttributeWarns)
{
    const char xml[] = "<val a=\"x\" b=\"y\"/>";
    auto ptree = Tests::readXml(xml);
    Callbacks cbs;
    {
        auto const conf = makeConfigTree(std::move(ptree), cbs);
        auto children = conf.getAllChildren();
        EXPECT_EQ(1, children.size());
        // Read only 'a'; 'b' is left unread.
        EXPECT_EQ("x",
                  children[0].second->getConfigAttribute<std::string>("a"));
        EXPECT_ERR_WARN(cbs, false, false);
    }
    // Attribute 'b' was not read -> warning.
    EXPECT_ERR_WARN(cbs, false, true);
}

// A content-less child (no immediate data, no attributes, no sub-children)
// left unread still warns: a returned child that is never dereferenced counts
// as unread regardless of content. This is what lets getAllChildren() flag a
// child whose tag name was mistyped that would otherwise be silently
// ignored.
TEST(BaseLibConfigTree, GetAllChildrenContentlessChildUnreadWarns)
{
    const char xml[] =
        "<a/>"
        "<b/>";
    auto ptree = Tests::readXml(xml);
    Callbacks cbs;
    {
        auto const conf = makeConfigTree(std::move(ptree), cbs);
        auto children = conf.getAllChildren();
        // Only inspect the tags; access neither child.
        EXPECT_EQ(2, children.size());
        EXPECT_EQ("a", children[0].first);
        EXPECT_EQ("b", children[1].first);
        EXPECT_ERR_WARN(cbs, false, false);
    }
    // Both content-less children went unread (never dereferenced) -> warning.
    EXPECT_ERR_WARN(cbs, false, true);
}

// Dereferencing a content-less child marks it read, so it is silent: touching
// the child is the signal that the tag was expected.
TEST(BaseLibConfigTree, GetAllChildrenContentlessChildAccessedStaysSilent)
{
    const char xml[] = "<a/>";
    auto ptree = Tests::readXml(xml);
    Callbacks cbs;
    {
        auto const conf = makeConfigTree(std::move(ptree), cbs);
        auto children = conf.getAllChildren();
        EXPECT_EQ(1, children.size());
        // Dereference the (content-less) child; nothing to read inside it.
        [[maybe_unused]] auto const& ct = *children[0].second;
        EXPECT_ERR_WARN(cbs, false, false);
    }
    // The child was accessed, so no unread-child warning.
    EXPECT_ERR_WARN(cbs, false, false);
}

// Partial consume, subtree flavour: descending into one child but leaving a
// sibling subtree untouched warns for the unread sibling, mirroring the scalar
// case in GetAllChildrenPartialConsumeWarnsForUnread.
TEST(BaseLibConfigTree, GetAllChildrenPartialConsumeSubtreeWarns)
{
    const char xml[] =
        "<outer1><inner>alpha</inner></outer1>"
        "<outer2><inner>beta</inner></outer2>";
    auto ptree = Tests::readXml(xml);
    Callbacks cbs;
    {
        auto const conf = makeConfigTree(std::move(ptree), cbs);
        auto children = conf.getAllChildren();
        EXPECT_EQ(2, children.size());
        // Descend into 'outer1' only; 'outer2' remains unread.
        auto inner1 = children[0].second->getConfigSubtree("inner");
        EXPECT_EQ("alpha", inner1.getValue<std::string>());
        EXPECT_ERR_WARN(cbs, false, false);
    }
    // 'outer2' (and its unread <inner>) -> warning.
    EXPECT_ERR_WARN(cbs, false, true);
}

// Lifetime: a Child owns its subtree (via the shared top-level pointer), so it
// may outlive the ConfigTree it came from. Here a single child is moved out and
// kept alive after the parent (and the Children vector) are destroyed; it is
// still readable and consuming it afterwards warns nothing.
TEST(BaseLibConfigTree, GetAllChildrenChildOutlivesParent)
{
    const char xml[] = "<a>1</a>";
    auto ptree = Tests::readXml(xml);
    Callbacks cbs;
    {
        std::optional<BaseLib::ConfigTree::Child> kept;
        {
            auto const conf = makeConfigTree(std::move(ptree), cbs);
            auto children = conf.getAllChildren();
            EXPECT_EQ(1, children.size());
            kept.emplace(std::move(children[0].second));
        }  // parent and Children vector destroyed here
        // Parent tore down cleanly: getAllChildren() marked the tag
        // consumed, and the moved-from entry left in the vector is inert.
        EXPECT_ERR_WARN(cbs, false, false);

        // The surviving child is still usable after its parent is gone.
        EXPECT_EQ("1", (*kept)->getValue<std::string>());
        EXPECT_ERR_WARN(cbs, false, false);
    }  // kept child destroyed here; its value was read -> no warning
    EXPECT_ERR_WARN(cbs, false, false);
}

TEST(BaseLibConfigTree, ChildLivesOnIfParentDies)
{
    const char xml[] =
        "<s>test</s>"
        "<t>Test</t>"
        "<u>data</u>";
    auto ptree = Tests::readXml(xml);

    Callbacks cbs;

    {
        std::optional<BaseLib::ConfigTree> opt_child;

        {
            auto const parent = makeConfigTree(std::move(ptree), cbs);
            opt_child = parent.getConfigSubtree("s");

            // do something with parent after subtree access
            EXPECT_EQ("Test", parent.getConfigParameter<std::string>("t"));
        }  // parent goes out of scope

        EXPECT_ERR_WARN(cbs, false,
                        true);  // warning because <u> has not been read.

        // do something with the child
        EXPECT_EQ("test", opt_child->getValue<std::string>());
    }

    EXPECT_ERR_WARN(cbs, false, false);
}
