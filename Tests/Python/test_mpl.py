import inspect
import math

import numpy as np
import ogs.mpl as mpl
import pytest


def get_public_classes(module):
    return [
        cls
        for name, cls in inspect.getmembers(module, inspect.isclass)
        if cls.__module__ == module.__name__
        and not name.startswith("_")
        and name not in {"VariableArray", "Variable"}
    ]


@pytest.mark.parametrize("cls", get_public_classes(mpl))
def test_class_and_members_have_docstrings(cls):
    assert cls.__doc__, f"Missing docstring for class: {cls.__name__}"
    assert cls.__init__.__doc__, f"Missing __init__ docstring in class: {cls.__name__}"

    for name, member in inspect.getmembers(cls):
        if name.startswith("_"):
            continue
        doc = getattr(member, "__doc__", None)
        assert doc, f"Missing docstring for member '{name}' in class '{cls.__name__}'"


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

    # Resetting of spatial position values
    with pytest.raises(TypeError):
        mpl.SpatialPosition().node_id = None
    pos = mpl.SpatialPosition()
    pos.node_id = 2
    assert pos.node_id == 2

    with pytest.raises(TypeError):
        mpl.SpatialPosition().element_id = None
    pos = mpl.SpatialPosition()
    pos.element_id = 5
    assert pos.element_id == 5

    with pytest.raises(RuntimeError, match=r"[bB]ad[ _]optional[ _]access"):
        mpl.SpatialPosition().coordinates = None
    with pytest.raises(RuntimeError, match=r"Unable to cast.*"):
        mpl.SpatialPosition().coordinates = []
    pos = mpl.SpatialPosition()
    pos.coordinates = [1.0, 2.0, 3.0]
    np.testing.assert_allclose(pos.coordinates, [1.0, 2.0, 3.0])


def test_constant_and_property_interface():
    va = mpl.VariableArray()
    p = mpl.Constant("C", 5)
    x = mpl.SpatialPosition(coords=[0, 1, 2])

    assert p.value(va, x, 0, 0) == 5
    assert p.value(va, va, x, 0, 0) == 5
    assert p.dValue(va, mpl.Variable.temperature, x, 0, 0) == 0
    assert p.dValue(va, va, mpl.Variable.temperature, x, 0, 0) == 0


def test_linear_property():
    va = mpl.VariableArray()
    x = mpl.SpatialPosition(coords=[0, 1, 2])

    p = mpl.Linear("linear", 1000, [(mpl.Variable.temperature, 273.15, -4e-4)])
    va.temperature = 373.15
    assert p.value(va, x, 0, 0) == 960
    assert p.dValue(va, mpl.Variable.temperature, x, 0, 0) == 1000 * (-4e-4)
    assert p.dValue(va, mpl.Variable.liquid_phase_pressure, x, 0, 0) == 0


def test_bilinear_property():
    va = mpl.VariableArray()
    x = mpl.SpatialPosition(coords=[0, 1, 2])

    p = mpl.Linear(
        "bilinear",
        1000,
        [
            (mpl.Variable.temperature, 273.15, -4e-4),
            (mpl.Variable.liquid_phase_pressure, 1e5, 4.65e-10),
        ],
    )
    va.temperature = 373.15
    va.liquid_phase_pressure = 1e6
    assert p.value(va, x, 0, 0) == pytest.approx(960.4185, rel=1e-12)
    assert p.dValue(va, mpl.Variable.temperature, x, 0, 0) == 1000 * (-4e-4)
    assert p.dValue(va, mpl.Variable.liquid_phase_pressure, x, 0, 0) == 1000 * 4.65e-10

    va.temperature = 273.15
    assert p.value(va, x, 0, 0) == pytest.approx(1000.4185, rel=1e-12)
    assert p.dValue(va, mpl.Variable.temperature, x, 0, 0) == 1000 * (-4e-4)
    assert p.dValue(va, mpl.Variable.liquid_phase_pressure, x, 0, 0) == 1000 * 4.65e-10
