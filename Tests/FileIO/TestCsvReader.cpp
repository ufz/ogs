/**
 * \file   TestCsvReader.cpp
 * \author Karsten Rink
 * \date   2015-04-09
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include <cstdio>

#include "gtest/gtest.h"

#include "BaseLib/BuildInfo.h"
#include "Applications/FileIO/CsvInterface.h"
#include "GeoLib/Point.h"

class CsvInterfaceTest : public ::testing::Test
{
public:
    CsvInterfaceTest()
        : _file_name(BaseLib::BuildInfo::tests_tmp_path+"test.csv")
    {
        std::ofstream out(_file_name);
        out << "id\tx\ty\tz\tname\tvalue1\tvalue_two\n";
        out << "0\t642015.538\t5724666.445\t391.759\ttest_a\t11.05303121\t436.913\t133\n";
        out << "1\t642015.49\t724667.426\t391.85\ttest_b\t51.65503659\n";
        out << "2\t642015.379\t5724668.424\t391.914\ttest_c\t437.068\t135\t2\n";
        out << "3\t642015.318\t5724669.411\t392.033\ttest_d\t51.65505447\t11.05302923\n";
        out << "4\t642015.275\t5724670.403\t392.172\ttest_e\t437.326\t137\t392.172\n";
        out << "5\t642015.288\t5724671.407\t392.232\ttest_f\n";
        out << "6\t642015.231\t5724672.403\t392.281\ttest_g\t\t437.435\n";
        out << "7\t642015.232\t5724673.384\t392.385\ttest_h\t11.05302961\t437.539\n";
        out << "8\t642015.153\t5724674.372\t392.428\ttest_i\t51.65509909\t11.05302887\n";
        out << "9\t642015.137\t5724675.377\t392.485\ttest_j\t51.65510812\t11.05302905\n";
        out.close();
    }

    ~CsvInterfaceTest() override { std::remove(_file_name.c_str()); }

protected:
    int _result;
    std::string _file_name;
};

/// Reading 3D points
TEST_F(CsvInterfaceTest, SimpleReadPoints)
{
    std::vector<GeoLib::Point*> points;
    std::vector<GeoLib::Point*> points2;
    _result = FileIO::CsvInterface::readPoints(_file_name, '\t', points);
    ASSERT_EQ(0, _result);
    _result = FileIO::CsvInterface::readPoints(_file_name, '\t', points2, "x", "y", "z");
    ASSERT_EQ(0, _result);
    ASSERT_TRUE(points.size() == 10);
    ASSERT_TRUE(points2.size() == 10);
    for (std::size_t i=0; i<points.size(); ++i)
    {
        ASSERT_TRUE((*points[i])[1] == (*points2[i])[0]);
        ASSERT_TRUE((*points[i])[2] == (*points2[i])[1]);
    }
    for (auto p : points)
        delete p;
    for (auto p : points2)
        delete p;
}

/// Dealing with unconvertable data types
TEST_F(CsvInterfaceTest, StringInPointColumn)
{
    std::vector<GeoLib::Point*> points;
    _result = FileIO::CsvInterface::readPoints(_file_name, '\t', points, "x", "y", "name");
    ASSERT_EQ(10, _result);
    ASSERT_TRUE(points.empty());
}

/// Dealing with not existing columns
TEST_F(CsvInterfaceTest, WrongColumnName)
{
    std::vector<GeoLib::Point*> points;
    _result = FileIO::CsvInterface::readPoints(_file_name, '\t', points, "x", "y", "wrong_column_name");
    ASSERT_EQ(-1, _result);
    ASSERT_TRUE(points.empty());

    _result = FileIO::CsvInterface::readPoints(_file_name, '\t', points, "wrong_column_name", "y", "id");
    ASSERT_EQ(-1, _result);
    ASSERT_TRUE(points.empty());
}

/// Dealing with missing values
TEST_F(CsvInterfaceTest, MissingValues)
{
    std::vector<GeoLib::Point*> points;
    _result = FileIO::CsvInterface::readPoints(_file_name, '\t', points, "z", "value1", "value_two");
    ASSERT_EQ(3, _result);
    ASSERT_EQ(7, points.size());
    ASSERT_NEAR(437.539, (*points[4])[2], std::numeric_limits<double>::epsilon());
    for (auto p : points)
        delete p;
}

/// Reading 2D points
TEST_F(CsvInterfaceTest, Points2D)
{
    std::vector<GeoLib::Point*> points;
    _result = FileIO::CsvInterface::readPoints(_file_name, '\t', points, "x", "y");
    ASSERT_EQ(0, _result);
    ASSERT_EQ(10, points.size());
    for (auto & point : points)
        ASSERT_NEAR(0, (*point)[2], std::numeric_limits<double>::epsilon());
    for (auto p : points)
        delete p;
}

/// Dealing with non-sequential column order
TEST_F(CsvInterfaceTest, CoordinateOrder)
{
    std::vector<GeoLib::Point*> points1;
    std::vector<GeoLib::Point*> points2;
    std::vector<GeoLib::Point*> points3;
    _result = FileIO::CsvInterface::readPoints(_file_name, '\t', points1, "id", "y", "z");
    ASSERT_EQ(0, _result);
    _result = FileIO::CsvInterface::readPoints(_file_name, '\t', points2, "id", "z", "y");
    ASSERT_EQ(0, _result);
    _result = FileIO::CsvInterface::readPoints(_file_name, '\t', points3, "y", "id", "z");
    ASSERT_EQ(0, _result);
    ASSERT_EQ(10, points1.size());
    ASSERT_EQ(10, points2.size());
    ASSERT_EQ(10, points3.size());
    for (std::size_t i=0; i<points1.size(); ++i)
    {
        ASSERT_EQ((*points1[i])[1], (*points2[i])[2]);
        ASSERT_EQ((*points1[i])[2], (*points2[i])[1]);
        ASSERT_EQ((*points3[i])[0], (*points2[i])[2]);
        ASSERT_EQ((*points3[i])[1], (*points2[i])[0]);
        ASSERT_EQ((*points1[i])[0], (*points3[i])[1]);
        ASSERT_EQ((*points1[i])[1], (*points3[i])[0]);
    }
    for (auto p : points1)
        delete p;
    for (auto p : points2)
        delete p;
    for (auto p : points3)
        delete p;
}

/// Getting single columns
TEST_F(CsvInterfaceTest, GetColumn)
{
    std::vector<std::string> names;
    _result = FileIO::CsvInterface::readColumn<std::string>(_file_name, '\t', names, "name");
    ASSERT_EQ(0, _result);
    ASSERT_EQ(10, names.size());

    std::vector<double> values;
    _result = FileIO::CsvInterface::readColumn<double>(_file_name, '\t', values, "value_two");
    ASSERT_EQ(2, _result);
    ASSERT_EQ(8, values.size());
}

/// Dealing with non-existing column
TEST_F(CsvInterfaceTest, NonExistingColumn)
{
    std::vector<double> values;
    _result = FileIO::CsvInterface::readColumn<double>(_file_name, '\t', values, "value2");
    ASSERT_EQ(-1, _result);
    ASSERT_TRUE(values.empty());
}

/// Dealing with wrong data type
TEST_F(CsvInterfaceTest, WrongDataType)
{
    std::vector<double> values;
    _result = FileIO::CsvInterface::readColumn<double>(_file_name, '\t', values, "name");
    ASSERT_EQ(10, _result);
    ASSERT_TRUE(values.empty());

    std::vector<std::string> names;
    _result = FileIO::CsvInterface::readColumn<std::string>(_file_name, '\t', names, "value1");
    ASSERT_EQ(2, _result);
    ASSERT_EQ(8, names.size());
}

