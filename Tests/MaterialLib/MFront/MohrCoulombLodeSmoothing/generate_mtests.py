# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

"""Generate Lode-angle MTests and print their names as a CMake list."""

import argparse
import re
from pathlib import Path


def material_property(name, value, kind="constant"):
    return f"@MaterialProperty<{kind}> '{name}' {value};"


def generate_mtests(library, output_directory):
    template = Path(__file__).with_name("MohrCoulombAbboSloan_lode.mtest.in")
    source = template.read_text(encoding="utf-8")
    transition_angle = 20
    # (name suffix, angle, lambda, radius factor): a plastic return inside the
    # transition and in the smoothing, and an elastic state just inside the
    # yield surface in the smoothing.
    cases = (
        ("", transition_angle - 0.01, "1e-5", "1"),
        ("", transition_angle + 5, "1e-5", "1"),
        ("-elastic", transition_angle + 5, "0", "1-2e-4"),
    )
    output_directory.mkdir(parents=True, exist_ok=True)
    test_names = []
    for variant in ("", "Aniso", "Ortho", "UBI", "UBIOrtho"):
        behaviour = f"MohrCoulombAbboSloan{variant}"
        elastic_properties = [
            material_property("YoungModulus", "'E'", "function"),
            material_property("PoissonRatio", "'nu'", "function"),
        ]
        if variant not in ("", "UBI"):
            elastic_properties = [
                material_property(f"YoungModulus{axis}", "'E'", "function")
                for axis in (1, 2, 3)
            ]
            for axes in (12, 23, 13):
                elastic_properties.extend(
                    (
                        material_property(f"PoissonRatio{axes}", "'nu'", "function"),
                        material_property(
                            f"ShearModulus{axes}", "'E/(2*(1+nu))'", "function"
                        ),
                    )
                )
            elastic_properties.append("")
        hill_properties = ""
        if variant == "Aniso":
            hill_properties = "\n".join(
                material_property(name, 1) for name in ("NormalFactor", "ShearFactor")
            )
        weak_plane_properties = ""
        plastic_strain_variable = "EquivalentPlasticStrain"
        if variant.startswith("UBI"):
            # Keep the weak plane elastic: this regression tests the matrix criterion.
            weak_plane_properties = "\n".join(
                (
                    material_property("PlaneCohesion", "1e6"),
                    material_property("PlaneFrictionAngle", "'phi_deg'", "function"),
                    material_property("PlaneDilatancyAngle", "'psi_deg'", "function"),
                )
            )
            plastic_strain_variable = "EquivalentPlasticStrainMatrix"
        for suffix, angle, plastic_multiplier, radius_factor in cases:
            for sign in (-1, 1):
                signed_angle = f"{sign * angle:g}"
                test_name = f"mtest-{behaviour}-lode{suffix}-{signed_angle}"
                values = {
                    "_ogs_mfront_behaviour_lib": library.as_posix(),
                    "_lode_transition_angle": str(transition_angle),
                    "_lode_lambda": plastic_multiplier,
                    "_lode_radius_factor": radius_factor,
                    "lode_behaviour": behaviour,
                    "_lode_elastic_properties": "\n".join(elastic_properties),
                    "_lode_hill_properties": hill_properties,
                    "_lode_weak_plane_properties": weak_plane_properties,
                    "_lode_plastic_strain_variable": plastic_strain_variable,
                    "_lode_sign": str(sign),
                    "_lode_angle": signed_angle,
                }
                rendered = re.sub(
                    r"@(\w+)@", lambda match, values=values: values[match[1]], source
                )
                (output_directory / f"{test_name}.mtest").write_text(
                    rendered, encoding="utf-8", newline="\n"
                )
                test_names.append(test_name)
    return test_names


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    try:
        names = generate_mtests(args.library, args.output_directory)
    except (OSError, KeyError) as error:
        parser.exit(1, f"Lode-angle MTest generation failed: {error}\n")
    print(";".join(names))


if __name__ == "__main__":
    main()
