/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include <cstdio>

#include "gtest/gtest.h"

#include <string>
#include <vector>

#include "BaseLib/BuildInfo.h"
#include "Applications/FileIO/CsvInterface.h"

TEST(CsvWriter, WriteReadTest)
{
    std::string test_file(BaseLib::BuildInfo::tests_tmp_path + "TestData.csv");

    std::vector<std::string> str_vec {"Red", "Orange", "Yellow", "Green", "Blue", "Indigo", "Violet" };
    std::vector<int> int_vec { 1, 2, 4, 8, 16, 32, 64 };
    std::vector<double> dbl_vec;
    std::srand ( static_cast<unsigned>(std::time(nullptr)) );
    for (std::size_t i=0; i<int_vec.size(); ++i)
        dbl_vec.push_back(static_cast<double>(std::rand()) / RAND_MAX);

    FileIO::CsvInterface csv;
    bool added;
    std::vector<std::string> vec_names { "String Vector", "Int Vector", "Double Vector"};
    added = csv.addVectorForWriting(vec_names[0], str_vec);
    ASSERT_TRUE(added);
    added = csv.addVectorForWriting(vec_names[1], int_vec);
    ASSERT_TRUE(added);
    added = csv.addVectorForWriting(vec_names[2], dbl_vec);
    ASSERT_TRUE(added);
    int_vec.push_back(128);
    added = csv.addVectorForWriting(vec_names[1], int_vec);
    ASSERT_FALSE(added);
    ASSERT_EQ(3, csv.getNArrays());
    csv.addIndexVectorForWriting(str_vec.size());
    ASSERT_EQ(4, csv.getNArrays());
    int result = csv.writeToFile(test_file);
    ASSERT_EQ(1, result);

    std::vector<std::string> str_result;
    result = FileIO::CsvInterface::readColumn<std::string>(test_file, '\t', str_result, vec_names[0]);
    ASSERT_EQ(0, result);
    ASSERT_EQ(str_vec.size(), str_result.size());

    std::vector<int> idx_result;
    result = FileIO::CsvInterface::readColumn<int>(test_file, '\t', idx_result, "Index");
    ASSERT_EQ(0, result);
    ASSERT_EQ(dbl_vec.size(), idx_result.size());

    std::vector<int> int_result;
    result = FileIO::CsvInterface::readColumn<int>(test_file, '\t', int_result, vec_names[1]);
    ASSERT_EQ(0, result);
    // testing for vector length -1 because it had increased previously when testing size requirements
    ASSERT_EQ(int_vec.size()-1, int_result.size());

    std::vector<double> dbl_result;
    result = FileIO::CsvInterface::readColumn<double>(test_file, '\t', dbl_result, vec_names[2]);
    ASSERT_EQ(0, result);
    ASSERT_EQ(dbl_vec.size(), dbl_result.size());

    for (std::size_t i=0; i<str_vec.size(); ++i)
    {
        ASSERT_EQ(str_vec[i], str_result[i]);
        ASSERT_EQ(int_vec[i], int_result[i]);
        ASSERT_NEAR(dbl_vec[i], dbl_result[i], 10 * std::numeric_limits<double>::epsilon());
        ASSERT_EQ(static_cast<int>(i), idx_result[i]);
    }

    std::remove(test_file.c_str());
}

