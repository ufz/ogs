# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from xml.dom import minidom


def add_parameter_property(properties, name_, parameter_name):
    new_property = ET.SubElement(properties, "property")
    ET.SubElement(new_property, "name").text = name_
    ET.SubElement(new_property, "type").text = "Parameter"
    ET.SubElement(new_property, "parameter_name").text = parameter_name


def add_liquid_phase(phases, parameter_name_viscosity, parameter_name_density):
    phase = ET.SubElement(phases, "phase")
    ET.SubElement(phase, "type").text = "AqueousLiquid"
    properties = ET.SubElement(phase, "properties")

    add_parameter_property(properties, "viscosity", parameter_name_viscosity)
    add_parameter_property(properties, "density", parameter_name_density)


def add_solid_phase(phases, parameter_name_density):
    phase = ET.SubElement(phases, "phase")
    ET.SubElement(phase, "type").text = "Solid"
    properties = ET.SubElement(phase, "properties")

    add_parameter_property(properties, "density", parameter_name_density)


def add_medium_properties_common(properties, storage_parameter, biot_parameter):
    add_parameter_property(properties, "storage", storage_parameter)
    add_parameter_property(properties, "biot_coefficient", biot_parameter)

    property_T = ET.SubElement(properties, "property")
    ET.SubElement(property_T, "name").text = "reference_temperature"
    ET.SubElement(property_T, "type").text = "Constant"
    ET.SubElement(property_T, "value").text = "293.15"


def add_medium_properties_matrix(
    medium,
    storage_parameter,
    biot_parameter,
    permeability_parameter,
    porosity_parameter,
):
    properties = ET.SubElement(medium, "properties")
    add_parameter_property(properties, "permeability", permeability_parameter)
    add_parameter_property(properties, "porosity", porosity_parameter)

    add_medium_properties_common(properties, storage_parameter, biot_parameter)


def add_medium_properties_fracture(
    medium, storage_parameter, biot_parameter, permeability_parameter
):
    properties = ET.SubElement(medium, "properties")
    property_kf = ET.SubElement(properties, "property")
    ET.SubElement(property_kf, "name").text = "permeability"
    if permeability_parameter == "0.0":
        ET.SubElement(property_kf, "type").text = "CubicLawPermeability"
    else:
        ET.SubElement(property_kf, "type").text = "Constant"
        ET.SubElement(property_kf, "value").text = permeability_parameter

    add_medium_properties_common(properties, storage_parameter, biot_parameter)


def add_matrix_medium(
    media,
    mat_id,
    k_m_tag,
    S_m_tag,
    biot_m_tag,
    phi_m_tag,
    mu_tag,
    rho_fr_tag,
    rho_sr_tag,
):
    medium = ET.SubElement(media, "medium", id=mat_id)
    phases = ET.SubElement(medium, "phases")

    add_liquid_phase(phases, mu_tag.text, rho_fr_tag.text)

    add_solid_phase(phases, rho_sr_tag.text)

    add_medium_properties_matrix(
        medium,
        storage_parameter=S_m_tag.text,
        biot_parameter=biot_m_tag.text,
        permeability_parameter=k_m_tag.text,
        porosity_parameter=phi_m_tag.text,
    )


def add_fracture_medium(
    media,
    mat_id,
    if_constant_permeability,
    fracture_props_tag,
    mu_tag,
    rho_fr_tag,
    rho_sr_tag,
):
    medium = ET.SubElement(media, "medium", id=mat_id)
    phases = ET.SubElement(medium, "phases")

    add_liquid_phase(phases, mu_tag.text, rho_fr_tag.text)
    add_solid_phase(phases, rho_sr_tag.text)

    S_f_tag = fracture_props_tag.find("specific_storage")
    fracture_props_tag.remove(S_f_tag)
    b_f_tag = fracture_props_tag.find("biot_coefficient")
    fracture_props_tag.remove(b_f_tag)
    k_f_tag = fracture_props_tag.find("permeability_model")
    fracture_props_tag.remove(k_f_tag)

    add_medium_properties_fracture(
        medium,
        storage_parameter=S_f_tag.text,
        biot_parameter=b_f_tag.text,
        permeability_parameter=if_constant_permeability,
    )


def prettify(elem):
    """Return a pretty-printed XML string for the Element without extra blank lines."""
    encoding = "ISO-8859-1"
    rough_string = ET.tostring(elem, encoding=encoding)
    reparsed = minidom.parseString(rough_string)
    pretty_xml = reparsed.toprettyxml(indent="    ", encoding=encoding)

    pretty_xml_decode = pretty_xml.decode(encoding)
    # Remove extra blank lines
    return "\n".join([line for line in pretty_xml_decode.splitlines() if line.strip()])


def insert_media_tags_after_time_loop(input_file, output_file):
    parser = ET.XMLParser(target=ET.TreeBuilder(insert_comments=True))
    tree = ET.parse(input_file, parser=parser)
    root = tree.getroot()

    if root.find(".//media") is not None:
        print(f"The project file {input_file} has already been processed")
        return False

    process_tag = root.find(".//processes/process")
    mu_tag = process_tag.find("fluid_viscosity")
    process_tag.remove(mu_tag)
    rho_fr_tag = process_tag.find("fluid_density")
    process_tag.remove(rho_fr_tag)
    rho_sr_tag = process_tag.find("solid_density")
    process_tag.remove(rho_sr_tag)

    k_m_tag = process_tag.find("intrinsic_permeability")
    process_tag.remove(k_m_tag)
    S_m_tag = process_tag.find("specific_storage")
    process_tag.remove(S_m_tag)
    biot_m_tag = process_tag.find("biot_coefficient")
    process_tag.remove(biot_m_tag)
    phi_m_tag = process_tag.find("porosity")
    process_tag.remove(phi_m_tag)

    constitutive_relation_tags = process_tag.findall("constitutive_relation")

    matrix_ids = []
    if len(constitutive_relation_tags) > 1:
        for tag in constitutive_relation_tags:
            matrix_ids.append(tag.attrib["id"])
    else:
        matrix_ids.append("0")

    fracture_props_tag = process_tag.find("fracture_properties")
    fracture_id_tag = fracture_props_tag.find("material_id")
    fracture_id = fracture_id_tag.text
    if len(matrix_ids) > 1:
        if fracture_id in matrix_ids:
            print(
                f"The fracture material id {fracture_id} is  the same as one"
                "of one matrix material IDs."
            )
    else:
        if int(fracture_id) == 0:
            matrix_ids[0] = "1"
            fracture_id = "0"
    fracture_props_tag.remove(fracture_id_tag)

    time_loop = root.find(".//time_loop")

    if time_loop is None:
        print("<time_loop> tag not found in the project file.")
        return False

    # Create the new <media> structure
    media = ET.Element("media")

    # Matrix medium
    for matrix_id in matrix_ids:
        add_matrix_medium(
            media,
            matrix_id,
            k_m_tag,
            S_m_tag,
            biot_m_tag,
            phi_m_tag,
            mu_tag,
            rho_fr_tag,
            rho_sr_tag,
        )

    # Fracture medium
    fracture_permeability_tag = fracture_props_tag.find("permeability_model")
    fracture_permeability_type = fracture_permeability_tag.find("type").text
    const_k_f = (
        fracture_permeability_tag.find("value").text
        if fracture_permeability_type == "ConstantPermeability"
        else "0.0"
    )

    add_fracture_medium(
        media,
        fracture_id,
        const_k_f,
        fracture_props_tag,
        mu_tag,
        rho_fr_tag,
        rho_sr_tag,
    )

    # Insert the <media> element after <time_loop>
    parent = root
    index = list(parent).index(time_loop)
    parent.insert(index + 1, media)

    with Path(output_file).open("w") as f:
        f.write(prettify(root))
        f.write("\n")

    return True


def main():
    if "--help" in sys.argv:
        print(
            "This script converts the LIE_HM project file to the new "
            "output syntax for MPL properties."
        )
        print("Usage: convert_LIE_HM2MPL.py <input_file> <output_file>")
        return

    if len(sys.argv) != 3:
        print("Usage: convert_LIE_HM2MPL <input_file> <output_file>")
        return

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    try:
        status = insert_media_tags_after_time_loop(input_file, output_file)
        if status:
            print("The conversion succeeds!")

    except FileNotFoundError:
        print(f"Error: The file {input_file} is not found")
    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    main()
