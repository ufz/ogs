# ruff: noqa: E402
import pytest

ogstools = pytest.importorskip("ogstools")

import shutil
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import ogstools as ot
from ogs import cli

refrigerant_density = 992.92  # kg m^-3
refrigerant_heat_capacity = 4068  # J kg^-1 K^-1


# [:,0] coords and [:,1] values
heat_power = np.array([[0, 300, 600], [-300, -400, -200]])
cop_heating = np.array([[268.15, 293.15], [4, 6.5]])
dhw_power = np.array([[0, 300, 600], [-100, -150, -250]])
cop_dhw = np.array([[268.15, 293.15], [2.5, 4.7]])
cool_power = np.array([[0, 300, 600], [30, 60, 120]])
cop_cool = np.array([[283.15, 293.15], [7.5, 6.1]])
flow_curve = np.array([[0, 300, 600], [2e-4, 2e-4, 2e-4]])

srcdir = Path(__file__).parent.parent.parent
testsrcdir = srcdir / "Tests/Data/Parabolic/T/3D_Beier_sandbox"

prj = "beier_sandbox_advanced_building_power_curves.prj"
base_prj_path = (
    testsrcdir / prj
)  # case buildingPowerCurveHotWaterCurveActiveCoolingCurveFlowCurve

bhe_top_point = (0, 2.5, 0)


def buildingPowerCurveHotWaterCurveActiveCoolingCurveFlowCurve(test_dir):
    def test(outflow_temperature, time):
        power_heating = np.interp(time, heat_power[0, :], heat_power[1, :])
        cop_heat = np.interp(outflow_temperature, cop_heating[0, :], cop_heating[1, :])
        power_hot_water = np.interp(time, dhw_power[0, :], dhw_power[1, :])
        cop_hot_water = np.interp(outflow_temperature, cop_dhw[0, :], cop_dhw[1, :])
        power_cooling = np.interp(time, cool_power[0, :], cool_power[1, :])
        cop_cooling = np.interp(outflow_temperature, cop_cool[0, :], cop_cool[1, :])
        flow_rate = np.interp(time, flow_curve[0, :], flow_curve[1, :])

        power = (
            power_heating
            - power_heating / cop_heat
            + power_hot_water
            - power_hot_water / cop_hot_water
            + power_cooling
            - power_cooling / cop_cooling
        )

        return power / (flow_rate * refrigerant_density * refrigerant_heat_capacity)

    prj_path = test_dir / prj
    # return the modified prj and a callable to check the result in dependence of temperature and time
    return prj_path, test


def buildingPowerCurveHotWaterCurvePassiveCoolingCurveFlowCurve(test_dir):
    def test(outflow_temperature, time):
        power_heating = np.interp(time, heat_power[0, :], heat_power[1, :])
        cop_heat = np.interp(outflow_temperature, cop_heating[0, :], cop_heating[1, :])
        power_hot_water = np.interp(time, dhw_power[0, :], dhw_power[1, :])
        cop_hot_water = np.interp(outflow_temperature, cop_dhw[0, :], cop_dhw[1, :])
        power_cooling = np.interp(time, cool_power[0, :], cool_power[1, :])
        flow_rate = np.interp(time, flow_curve[0, :], flow_curve[1, :])

        power = (
            power_heating
            - power_heating / cop_heat
            + power_hot_water
            - power_hot_water / cop_hot_water
            + power_cooling
        )

        return power / (flow_rate * refrigerant_density * refrigerant_heat_capacity)

    prj_path = (
        test_dir / "buildingPowerCurveHotWaterCurvePassiveCoolingCurveFlowCurve.prj"
    )

    model = ot.Project(input_file=base_prj_path, output_file=prj_path)

    model.replace_text(
        value="false",
        xpath="./processes/process/borehole_heat_exchangers/borehole_heat_exchanger/flow_and_temperature_control/cooling/active",
    )

    model.remove_element(
        xpath="./processes/process/borehole_heat_exchangers/borehole_heat_exchanger/flow_and_temperature_control/cooling/cop_curve",
    )
    model.write_input()
    return prj_path, test


def buildingPowerCurveHotWaterCurveFlowCurve(test_dir):
    def test(outflow_temperature, time):
        power_heating = np.interp(time, heat_power[0, :], heat_power[1, :])
        cop_heat = np.interp(outflow_temperature, cop_heating[0, :], cop_heating[1, :])
        power_hot_water = np.interp(time, dhw_power[0, :], dhw_power[1, :])
        cop_hot_water = np.interp(outflow_temperature, cop_dhw[0, :], cop_dhw[1, :])
        flow_rate = np.interp(time, flow_curve[0, :], flow_curve[1, :])

        power = (
            power_heating
            - power_heating / cop_heat
            + power_hot_water
            - power_hot_water / cop_hot_water
        )

        return power / (flow_rate * refrigerant_density * refrigerant_heat_capacity)

    prj_path = test_dir / "buildingPowerCurveHotWaterCurveFlowCurve.prj"

    model = ot.Project(input_file=base_prj_path, output_file=prj_path)

    model.remove_element(
        xpath="./processes/process/borehole_heat_exchangers/borehole_heat_exchanger/flow_and_temperature_control/cooling",
    )
    model.write_input()
    return prj_path, test


def buildingPowerCurveActiveCoolingCurveFlowCurve(test_dir):
    def test(outflow_temperature, time):
        power_heating = np.interp(time, heat_power[0, :], heat_power[1, :])
        cop_heat = np.interp(outflow_temperature, cop_heating[0, :], cop_heating[1, :])
        power_heating = np.interp(time, heat_power[0, :], heat_power[1, :])
        power_cooling = np.interp(time, cool_power[0, :], cool_power[1, :])
        cop_cooling = np.interp(outflow_temperature, cop_cool[0, :], cop_cool[1, :])
        flow_rate = np.interp(time, flow_curve[0, :], flow_curve[1, :])

        power = (
            power_heating
            - power_heating / cop_heat
            + power_cooling
            - power_cooling / cop_cooling
        )

        return power / (flow_rate * refrigerant_density * refrigerant_heat_capacity)

    prj_path = test_dir / "buildingPowerCurveActiveCoolingCurveFlowCurve.prj"

    model = ot.Project(input_file=base_prj_path, output_file=prj_path)

    model.remove_element(
        xpath="./processes/process/borehole_heat_exchangers/borehole_heat_exchanger/flow_and_temperature_control/hot_water",
    )
    model.write_input()
    return prj_path, test


def buildingPowerCurvePassiveCoolingCurveFlowCurve(test_dir):
    def test(outflow_temperature, time):
        power_heating = np.interp(time, heat_power[0, :], heat_power[1, :])
        cop_heat = np.interp(outflow_temperature, cop_heating[0, :], cop_heating[1, :])
        power_cooling = np.interp(time, cool_power[0, :], cool_power[1, :])
        flow_rate = np.interp(time, flow_curve[0, :], flow_curve[1, :])

        power = power_heating - power_heating / cop_heat + power_cooling

        return power / (flow_rate * refrigerant_density * refrigerant_heat_capacity)

    prj_path = test_dir / "buildingPowerCurvePassiveCoolingCurveFlowCurve.prj"

    model = ot.Project(input_file=base_prj_path, output_file=prj_path)

    model.remove_element(
        xpath="./processes/process/borehole_heat_exchangers/borehole_heat_exchanger/flow_and_temperature_control/hot_water",
    )
    model.replace_text(
        value="false",
        xpath="./processes/process/borehole_heat_exchangers/borehole_heat_exchanger/flow_and_temperature_control/cooling/active",
    )

    model.remove_element(
        xpath="./processes/process/borehole_heat_exchangers/borehole_heat_exchanger/flow_and_temperature_control/cooling/cop_curve",
    )

    model.write_input()
    return prj_path, test


def hotWaterCurveActiveCoolingCurveFlowCurve(test_dir):
    def test(outflow_temperature, time):
        power_hot_water = np.interp(time, dhw_power[0, :], dhw_power[1, :])
        cop_hot_water = np.interp(outflow_temperature, cop_dhw[0, :], cop_dhw[1, :])
        power_cooling = np.interp(time, cool_power[0, :], cool_power[1, :])
        cop_cooling = np.interp(outflow_temperature, cop_cool[0, :], cop_cool[1, :])
        flow_rate = np.interp(time, flow_curve[0, :], flow_curve[1, :])

        power = (
            power_hot_water
            - power_hot_water / cop_hot_water
            + power_cooling
            - power_cooling / cop_cooling
        )

        return power / (flow_rate * refrigerant_density * refrigerant_heat_capacity)

    prj_path = test_dir / "hotWaterCurveActiveCoolingCurveFlowCurve.prj"

    model = ot.Project(input_file=base_prj_path, output_file=prj_path)

    model.remove_element(
        xpath="./processes/process/borehole_heat_exchangers/borehole_heat_exchanger/flow_and_temperature_control/heating",
    )
    model.write_input()
    return prj_path, test


def hotWaterCurvePassiveCoolingCurveFlowCurve(test_dir):
    def test(outflow_temperature, time):
        power_hot_water = np.interp(time, dhw_power[0, :], dhw_power[1, :])
        cop_hot_water = np.interp(outflow_temperature, cop_dhw[0, :], cop_dhw[1, :])
        power_cooling = np.interp(time, cool_power[0, :], cool_power[1, :])
        flow_rate = np.interp(time, flow_curve[0, :], flow_curve[1, :])

        power = power_hot_water - power_hot_water / cop_hot_water + power_cooling

        return power / (flow_rate * refrigerant_density * refrigerant_heat_capacity)

    prj_path = test_dir / "hotWaterCurvePassiveCoolingCurveFlowCurve.prj"

    model = ot.Project(input_file=base_prj_path, output_file=prj_path)

    model.remove_element(
        xpath="./processes/process/borehole_heat_exchangers/borehole_heat_exchanger/flow_and_temperature_control/heating",
    )
    model.replace_text(
        value="false",
        xpath="./processes/process/borehole_heat_exchangers/borehole_heat_exchanger/flow_and_temperature_control/cooling/active",
    )

    model.remove_element(
        xpath="./processes/process/borehole_heat_exchangers/borehole_heat_exchanger/flow_and_temperature_control/cooling/cop_curve",
    )

    model.write_input()
    return prj_path, test


def buildingPowerCurveFlowCurve(test_dir):
    def test(outflow_temperature, time):
        power_heating = np.interp(time, heat_power[0, :], heat_power[1, :])
        cop_heat = np.interp(outflow_temperature, cop_heating[0, :], cop_heating[1, :])
        flow_rate = np.interp(time, flow_curve[0, :], flow_curve[1, :])

        power = power_heating - power_heating / cop_heat

        return power / (flow_rate * refrigerant_density * refrigerant_heat_capacity)

    prj_path = test_dir / "buildingPowerCurveFlowCurve.prj"

    model = ot.Project(input_file=base_prj_path, output_file=prj_path)

    model.remove_element(
        xpath="./processes/process/borehole_heat_exchangers/borehole_heat_exchanger/flow_and_temperature_control/hot_water",
    )
    model.remove_element(
        xpath="./processes/process/borehole_heat_exchangers/borehole_heat_exchanger/flow_and_temperature_control/cooling",
    )
    model.write_input()
    return prj_path, test


def hotWaterCurveFlowCurve(test_dir):
    def test(outflow_temperature, time):
        power_hot_water = np.interp(time, dhw_power[0, :], dhw_power[1, :])
        cop_hot_water = np.interp(outflow_temperature, cop_dhw[0, :], cop_dhw[1, :])
        flow_rate = np.interp(time, flow_curve[0, :], flow_curve[1, :])

        power = power_hot_water - power_hot_water / cop_hot_water

        return power / (flow_rate * refrigerant_density * refrigerant_heat_capacity)

    prj_path = test_dir / "hotWaterCurveFlowCurve.prj"

    model = ot.Project(input_file=base_prj_path, output_file=prj_path)

    model.remove_element(
        xpath="./processes/process/borehole_heat_exchangers/borehole_heat_exchanger/flow_and_temperature_control/heating",
    )
    model.remove_element(
        xpath="./processes/process/borehole_heat_exchangers/borehole_heat_exchanger/flow_and_temperature_control/cooling",
    )
    model.write_input()
    return prj_path, test


def activeCoolingCurveFlowCurve(test_dir):
    def test(outflow_temperature, time):
        power_cooling = np.interp(time, cool_power[0, :], cool_power[1, :])
        cop_cooling = np.interp(outflow_temperature, cop_cool[0, :], cop_cool[1, :])
        flow_rate = np.interp(time, flow_curve[0, :], flow_curve[1, :])

        power = power_cooling - power_cooling / cop_cooling

        return power / (flow_rate * refrigerant_density * refrigerant_heat_capacity)

    prj_path = test_dir / "activeCoolingCurveFlowCurve.prj"

    model = ot.Project(input_file=base_prj_path, output_file=prj_path)

    model.remove_element(
        xpath="./processes/process/borehole_heat_exchangers/borehole_heat_exchanger/flow_and_temperature_control/heating",
    )
    model.remove_element(
        xpath="./processes/process/borehole_heat_exchangers/borehole_heat_exchanger/flow_and_temperature_control/hot_water",
    )

    model.write_input()
    return prj_path, test


def passiveCoolingCurveFlowCurve(test_dir):
    def test(_, time):
        power_cooling = np.interp(time, cool_power[0, :], cool_power[1, :])
        flow_rate = np.interp(time, flow_curve[0, :], flow_curve[1, :])

        return power_cooling / (
            flow_rate * refrigerant_density * refrigerant_heat_capacity
        )

    prj_path = test_dir / "passiveCoolingCurveFlowCurve.prj"

    model = ot.Project(input_file=base_prj_path, output_file=prj_path)

    model.remove_element(
        xpath="./processes/process/borehole_heat_exchangers/borehole_heat_exchanger/flow_and_temperature_control/heating",
    )
    model.remove_element(
        xpath="./processes/process/borehole_heat_exchangers/borehole_heat_exchanger/flow_and_temperature_control/hot_water",
    )
    model.replace_text(
        value="false",
        xpath="./processes/process/borehole_heat_exchangers/borehole_heat_exchanger/flow_and_temperature_control/cooling/active",
    )

    model.remove_element(
        xpath="./processes/process/borehole_heat_exchangers/borehole_heat_exchanger/flow_and_temperature_control/cooling/cop_curve",
    )

    model.write_input()
    return prj_path, test


@pytest.fixture(scope="module")
def temp_test_dir():
    with TemporaryDirectory() as tmp_dir:
        test_dir = shutil.copytree(testsrcdir, Path(tmp_dir) / testsrcdir.name)
        yield test_dir


@pytest.mark.ogs_needs_serial_build
@pytest.mark.parametrize(
    "case_gen",
    [
        buildingPowerCurveHotWaterCurveActiveCoolingCurveFlowCurve,
        buildingPowerCurveHotWaterCurvePassiveCoolingCurveFlowCurve,
        buildingPowerCurveHotWaterCurveFlowCurve,
        buildingPowerCurveActiveCoolingCurveFlowCurve,
        buildingPowerCurvePassiveCoolingCurveFlowCurve,
        hotWaterCurveActiveCoolingCurveFlowCurve,
        hotWaterCurvePassiveCoolingCurveFlowCurve,
        buildingPowerCurveFlowCurve,
        hotWaterCurveFlowCurve,
        activeCoolingCurveFlowCurve,
        passiveCoolingCurveFlowCurve,
    ],
)
def test_bhe_advanced_building_power_curves(temp_test_dir, case_gen):
    prj_path, check_results = case_gen(temp_test_dir)

    status = cli.ogs(str(prj_path), o=temp_test_dir)
    assert status == 0  # OGS run successful

    results = ot.MeshSeries(temp_test_dir / "beier_sandbox.pvd")

    bhe_vector = ot.variables.temperature_BHE
    bhe_vector.output_unit = "K"

    bhe_temperatures = results.probe(
        bhe_top_point, bhe_vector[1, ["in", "out"]], interp_method="nearest"
    )
    timevalues = results.timevalues

    expected_bhe_temperature_difference = check_results(
        bhe_temperatures[1:, 1], timevalues[1:]
    )

    calculated_bhe_temperature_difference = bhe_temperatures[1:, 0] - (
        bhe_temperatures[1:, 1]
    )

    assert np.allclose(
        expected_bhe_temperature_difference,
        calculated_bhe_temperature_difference,
    )
