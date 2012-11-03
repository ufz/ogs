/**
 * \file
 * \author Dmitrij Naumov
 * \date   Nov. 2012
 * \brief  Testing the specialized version of quicksort.
 *
 * \copyright
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "gtest.h"
#include "quickcheck/quickcheck.hh"
#include "quicksort.h"

#include <algorithm>
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

