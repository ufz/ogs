/**
 * \file
 * \author Dmitrij Naumov
 * \date   Nov. 2012
 * \brief  Testing the specialized version of quicksort.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "gtest/gtest.h"
#include "quickcheck/quickcheck.hh"
#include "quicksort.h"

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

using namespace quickcheck;

// Quicksort result is sorted.
class QuicksortSortsAsSTLSort : public Property<std::vector<int>> {
    bool holdsFor(const std::vector<int>& xs)
    {
        std::vector<size_t> perm(xs.size());
        BaseLib::quicksort((int*)&(xs[0]), 0, xs.size(), &(perm[0]));
        return std::is_sorted(xs.begin(), xs.end());
    }

	bool isTrivial(const std::vector<int>& xs)
	{
		return xs.size() < 2;
	}

	const std::string classify(const std::vector<int>& xs)
	{
		std::string classification;
		if (xs.size() < 4)
			classification += "short";
		else
			classification += "long";

		if (std::is_sorted(xs.begin(), xs.end()))
			classification += ", sorted";
		else
			classification += ", unsorted";

		return classification;
	}
};

TEST(BaseLib, QuicksortSortsCorrectly) {
    QuicksortSortsAsSTLSort quicksortSortsAsSTLSort;
    ASSERT_TRUE(quicksortSortsAsSTLSort.check());
}

// Permutations of sorted, unique vector remain untouched.
class QuicksortCheckPermutations : public Property<std::vector<int>> {
	bool accepts(const std::vector<int>& xs)
	{
		return xs.size() > 2 &&
			std::is_sorted(xs.begin(), xs.end()) &&
			(std::adjacent_find(xs.begin(), xs.end()) == xs.end());
	}

    bool holdsFor(const std::vector<int>& xs)
    {
        std::vector<size_t> perm(xs.size());
		for (size_t i = 0; i < perm.size(); ++i)
			perm[i] = i;

        BaseLib::quicksort((int*)&(xs[0]), 0, xs.size(), &(perm[0]));

		for (size_t i = 0; i < perm.size(); ++i)
			if (perm[i] != i)
				return false;
        return true;
    }

	const std::string classify(const std::vector<int>& xs)
	{
		std::stringstream ss;
		ss << "size " << xs.size();
		return ss.str();
	}
};

TEST(BaseLib, QuicksortReportCorrectPermutations) {
    QuicksortCheckPermutations quicksortCheckPermutations;
    ASSERT_TRUE(quicksortCheckPermutations.check(200, 100000));
}

// Permutations of reverse sorted, unique vector is also reversed.
class QuicksortCheckPermutationsReverse : public Property<std::vector<int>> {
	bool accepts(const std::vector<int>& xs)
	{
		return xs.size() > 2 &&
			std::is_sorted(xs.rbegin(), xs.rend()) &&
			(std::adjacent_find(xs.rbegin(), xs.rend()) == xs.rend());
	}

    bool holdsFor(const std::vector<int>& xs)
    {
        std::vector<size_t> perm(xs.size());
		for (size_t i = 0; i < perm.size(); ++i)
			perm[i] = i;

        BaseLib::quicksort((int*)&(xs[0]), 0, xs.size(), &(perm[0]));

		for (size_t i = 0; i < perm.size(); ++i)
			if (perm[i] != perm.size() - i - 1)
				return false;
        return true;
    }

	const std::string classify(const std::vector<int>& xs)
	{
		std::stringstream ss;
		ss << "size " << xs.size();
		return ss.str();
	}
};

TEST(BaseLib, QuicksortReportCorrectPermutationsReverse) {
    QuicksortCheckPermutationsReverse quicksortCheckPermutationsReverse;
    ASSERT_TRUE(quicksortCheckPermutationsReverse.check(200, 100000));
}
