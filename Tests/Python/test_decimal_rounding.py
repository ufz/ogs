import random

from ogs.dev.ogs_log_summary import _round_up_2_digits as r


def test_without_change_upon_rounding():
    d_in = 0.3
    d_out = r(d_in)
    assert d_out == d_in

    d_in = 0.99
    d_out = r(d_in)
    assert d_out == d_in

    d_in = 2.3
    d_out = r(d_in)
    assert d_out == d_in

    d_in = 120
    d_out = r(d_in)
    assert d_out == d_in

    d_in = 0
    d_out = r(d_in)
    assert d_out == d_in

    d_in = -0.4
    d_out = r(d_in)
    assert d_out == -0.4

    d_in = -4.7
    d_out = r(d_in)
    assert d_out == d_in

    d_in = -27
    d_out = r(d_in)
    assert d_out == d_in

    d_in = -1300
    d_out = r(d_in)
    assert d_out == d_in


def test_with_change_upon_rounding():
    d_in = 0.4
    d_out = r(d_in)
    print(f"in:  {d_in:.17}")
    print(f"out: {d_out:.17}")
    # the above prints:
    # in:  0.40000000000000002
    # out: 0.40999999999999998
    # therefore the following holds:
    assert d_out == 0.41
    # Note: this corner case (and others) could probably be avoided by avoiding
    # the float to Decimal (i.e. binary to decimal) conversion altogether.

    d_in = 0.3999999
    d_out = r(d_in)
    assert d_out == 0.4

    d_in = 1.3999999
    d_out = r(d_in)
    assert d_out == 1.4

    d_in = 13.00001
    d_out = r(d_in)
    assert d_out == 14

    d_in = 13.99999
    d_out = r(d_in)
    assert d_out == 14

    d_in = -2.3
    d_out = r(d_in)
    print(f"in:  {d_in:.17}")
    print(f"out: {d_out:.17}")
    # the above prints:
    # in:  -2.2999999999999998
    # out: -2.2000000000000002
    # therefore the following holds:
    assert d_out == -2.2

    d_in = -0.41
    d_out = r(d_in)
    print(f"in:  {d_in:.17}")
    print(f"out: {d_out:.17}")
    # the above prints:
    # in:  -0.40999999999999998
    # out: -0.40000000000000002
    # therefore the following holds:
    assert d_out == -0.4

    d_in = -1.31
    d_out = r(d_in)
    assert d_out == -1.3

    d_in = -13.00001
    d_out = r(d_in)
    assert d_out == -13


def test_random():
    d_max = 1e9
    d_min = -d_max

    for _i in range(1000):
        d = random.uniform(d_min, d_max)
        d_up = r(d)

        assert d_up >= d  # rounding up
        assert (d_up >= 0) == (d >= 0)  # same sign
        # rounding "tightly" to two digits
        # note: 1 + eps will be rounded to 1.1, so 0.1 is the strict upper bound
        # for the relative change, since (1.1 - 1.0) / 1.0 == 0.1
        if d != 0:
            assert abs((d_up - d) / d) < 0.1
