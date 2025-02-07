import math

import numpy as np
import ogs.mpl as mpl
import pytest


def test_variable_array():
    va = mpl.VariableArray()
    assert math.isnan(va.capillary_pressure)
    va.capillary_pressure = 1e6
    assert va.capillary_pressure == 1e6

    assert va.deformation_gradient is None
    va.deformation_gradient = [0, 1, 2, 3, 4]
    assert np.allclose(va.deformation_gradient, [0, 1, 2, 3, 4])
    va.deformation_gradient = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    assert np.allclose(va.deformation_gradient, [0, 1, 2, 3, 4, 5, 6, 7, 8])

    assert va.mechanical_strain is None
    va.mechanical_strain = [0, 1, 2, 3]
    assert np.allclose(va.mechanical_strain, [0, 1, 2, 3])
    va.mechanical_strain = [0, 1, 2, 3, 4, 5]
    assert np.allclose(va.mechanical_strain, [0, 1, 2, 3, 4, 5])

    with pytest.raises(
        TypeError,
        match=r".*incompatible function arguments.*",
    ):
        va.capillary_pressure = [0, 1]
    with pytest.raises(
        TypeError,
        match=r".*incompatible function arguments.*",
    ):
        va.mechanical_strain = [0, 1, 2, 3, 4]
    with pytest.raises(
        TypeError,
        match=r".*incompatible function arguments.*",
    ):
        va.deformation_gradient = [0, 1, 2, 3]


def test_spatial_position():
    # Different coordinates constructions
    assert np.allclose(mpl.SpatialPosition(coords=[0, 1, 2]).coordinates, [0, 1, 2])
    assert np.allclose(mpl.SpatialPosition(coords=(0, 1, 2)).coordinates, [0, 1, 2])
    assert np.allclose(
        mpl.SpatialPosition(None, None, [0, 1, 2]).coordinates, [0, 1, 2]
    )
    assert np.allclose(
        mpl.SpatialPosition(None, None, (0, 1, 2)).coordinates, [0, 1, 2]
    )

    assert str(mpl.SpatialPosition(1, 2)) == "<SpatialPosition(1, 2, None)>"
    assert (
        str(mpl.SpatialPosition(coords=[0, 1, 2]))
        == "<SpatialPosition(None, None, [0. 1. 2.])>"
    )


def test_constant_and_property_interface():
    va = mpl.VariableArray()
    p = mpl.Constant("C", 5)
    x = mpl.SpatialPosition(coords=[0, 1, 2])

    assert p.value(va, x, 0, 0) == 5
    assert p.value(va, va, x, 0, 0) == 5
    assert p.dValue(va, mpl.Variable.temperature, x, 0, 0) == 0
    assert p.dValue(va, va, mpl.Variable.temperature, x, 0, 0) == 0


