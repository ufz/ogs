import ogstools.ogs6py as ogs6py


def addSolidPropertiesForHT(omit=""):
    if omit != "storage":
        model.media.add_property(
            medium_id="0",
            phase_type="Solid",
            name="storage",
            type="Constant",
            value="0.0",
        )
    if omit != "density":
        model.media.add_property(
            medium_id="0",
            phase_type="Solid",
            name="density",
            type="Constant",
            value="0.0",
        )
    if omit != "thermal_conductivity":
        model.media.add_property(
            medium_id="0",
            phase_type="Solid",
            name="thermal_conductivity",
            type="Constant",
            value="3.0",
        )
    if omit != "specific_heat_capacity":
        model.media.add_property(
            medium_id="0",
            phase_type="Solid",
            name="specific_heat_capacity",
            type="Constant",
            value="0.0",
        )


def addAqueousLiquidPropertiesForHT(omit=""):
    if omit != "specific_heat_capacity":
        model.media.add_property(
            medium_id="0",
            phase_type="AqueousLiquid",
            name="specific_heat_capacity",
            type="Constant",
            value="0.0",
        )
    if omit != "thermal_conductivity":
        model.media.add_property(
            medium_id="0",
            phase_type="AqueousLiquid",
            name="thermal_conductivity",
            type="Constant",
            value="0.65",
        )
    if omit != "density":
        model.media.add_property(
            medium_id="0",
            phase_type="AqueousLiquid",
            name="density",
            type="Constant",
            value="1e-3",
        )
    if omit != "viscosity":
        model.media.add_property(
            medium_id="0",
            phase_type="AqueousLiquid",
            name="viscosity",
            type="Constant",
            value="1e-3",
        )


def addMediumPropertiesForHT(omit=""):
    if omit != "thermal_longitudinal_dispersivity":
        model.media.add_property(
            medium_id="0",
            name="thermal_longitudinal_dispersivity",
            type="Constant",
            value="0.0",
        )
    if omit != "thermal_transversal_dispersivity":
        model.media.add_property(
            medium_id="0",
            name="thermal_transversal_dispersivity",
            type="Constant",
            value="0.0",
        )
    if omit != "permeability":
        model.media.add_property(
            medium_id="0", name="permeability", type="Constant", value="1e-14 0 0 1e-14"
        )
    if omit != "porosity":
        model.media.add_property(
            medium_id="0", name="porosity", type="Constant", value="1e-3"
        )


solid_properties = [
    "",
    "storage",
    "density",
    "thermal_conductivity",
    "specific_heat_capacity",
]
aqueousfluid_properties = [
    "",
    "specific_heat_capacity",
    "thermal_conductivity",
    "density",
    "viscosity",
]
medium_properties = [
    "",
    "thermal_longitudinal_dispersivity",
    "thermal_transversal_dispersivity",
    "permeability",
    "porosity",
]

for solid_property in solid_properties:
    for aqueousfluid_property in aqueousfluid_properties:
        for medium_property in medium_properties:
            if (
                solid_property == ""
                and aqueousfluid_property == ""
                and medium_property == ""
            ):
                continue
            model = ogs6py.Project(
                output_file="HT_"
                + solid_property
                + "_"
                + aqueousfluid_property
                + "_"
                + medium_property
                + ".prj"
            )
            model.mesh.add_mesh(filename="square_1x1_quad_1e3.vtu")
            model.geometry.add_geometry(filename="square_1x1.gml")
            model.processes.set_process(
                name="HT", type="HT", integration_order="2", specific_body_force="0 0"
            )
            model.processes.add_process_variable(
                process_variable="temperature", process_variable_name="temperature"
            )
            model.processes.add_process_variable(
                process_variable="pressure", process_variable_name="pressure"
            )
            model.processes.add_secondary_variable(
                internal_name="darcy_velocity", output_name="darcy_velocity"
            )

            addAqueousLiquidPropertiesForHT(aqueousfluid_property)
            addSolidPropertiesForHT(solid_property)
            addMediumPropertiesForHT(medium_property)

            model.time_loop.add_process(
                process="HT",
                nonlinear_solver_name="basic_picard",
                convergence_type="DeltaX",
                norm_type="NORM2",
                abstol="1e-15",
                time_discretization="BackwardEuler",
            )
            model.time_loop.set_stepping(
                process="HT",
                type="FixedTimeStepping",
                t_initial="0",
                t_end="1",
                repeat="4",
                delta_t="0.25",
            )
            model.time_loop.add_output(
                type="VTK",
                prefix="HT_test_",
                repeat="1",
                each_steps="10",
                variables=["temperature", "pressure"],
            )
            model.parameters.add_parameter(name="T0", type="Constant", value="0")
            model.parameters.add_parameter(name="P0", type="Constant", value="0")
            model.parameters.add_parameter(
                name="p_Dirichlet_left", type="Constant", value="1"
            )
            # model.parameters.add_parameter(name="p_Dirichlet_right", type="Constant", value="-1")
            model.parameters.add_parameter(
                name="t_Dirichlet_bottom", type="Constant", value="2"
            )
            # model.parameters.add_parameter(name="t_Dirichlet_top", type="Constant", value="1")
            model.process_variables.set_ic(
                process_variable_name="temperature",
                components="1",
                order="1",
                initial_condition="T0",
            )
            model.process_variables.add_bc(
                process_variable_name="temperature",
                geometrical_set="square_1x1_geometry",
                geometry="bottom",
                type="Dirichlet",
                component="0",
                parameter="t_Dirichlet_bottom",
            )
            # model.process_variables.add_bc(process_variable_name="temperature",
            #                        geometrical_set="square_1x1_geometry",
            #                        geometry="top",
            #                        type="Dirichlet",
            #                        component="1",
            #                        parameter="t_Dirichlet_top")
            model.process_variables.set_ic(
                process_variable_name="pressure",
                components="1",
                order="1",
                initial_condition="P0",
            )
            model.process_variables.add_bc(
                process_variable_name="pressure",
                geometrical_set="square_1x1_geometry",
                geometry="left",
                type="Dirichlet",
                component="1",
                parameter="p_Dirichlet_left",
            )
            # model.process_variables.addadd_bcBC(process_variable_name="pressure",
            #                        geometrical_set="square_1x1_geometry",
            #                        geometry="right",
            #                        type="Dirichlet",
            #                        component="1",
            #                        parameter="p_Dirichlet_right")
            model.nonlinear_solvers.add_non_lin_solver(
                name="basic_picard",
                type="Picard",
                max_iter="4",
                linear_solver="general_linear_solver",
            )
            model.linear_solvers.add_lin_solver(
                name="general_linear_solver",
                kind="lis",
                solver_type="cg",
                precon_type="jacobi",
                max_iteration_step="10000",
                error_tolerance="1e-16",
            )
            model.linear_solvers.add_lin_solver(
                name="general_linear_solver",
                kind="eigen",
                solver_type="CG",
                precon_type="DIAGONAL",
                max_iteration_step="10000",
                error_tolerance="1e-16",
            )
            model.linear_solvers.add_lin_solver(
                name="general_linear_solver",
                kind="petsc",
                solver_type="cg",
                precon_type="bjacobi",
                max_iteration_step="10000",
                error_tolerance="1e-16",
                prefix="foo",
            )
            model.write_input()
