"""
Script to search for OpenGeoSys PRJ files and move the
thermal_osmosis_coefficient property block from the solid phase to the medium
level, optionally converting it to thermal_osmosis_permeability.

Usage: python move_thermal_osmosis.py [ROOT_DIR] [--dry-run] [--backup] [--convert]

Options:
  --dry-run   Show what would be changed without making modifications
  --backup    Create backup of modified files
  --convert   Rename the property to thermal_osmosis_permeability and apply the
              conversion epsilon_T = kT * (mu / k) before writing

Scalar and diagonal tensor values are supported. Without --convert the property
element is moved verbatim, so any value survives unchanged. With --convert the
conversion is applied per diagonal component, which is only equivalent to the
relation kT = epsilon_T k / mu for diagonal kT and k; non-diagonal values are
rejected instead of silently mis-converted. epsilon_T is a scalar, so a
conversion whose components come out different is rejected as well.

The written file is not re-indented; run
`xmlstarlet format -s 4 old.prj > new.prj` afterwards.

Exits 1 if any file was refused or could not be parsed, and 0 otherwise,
including when no file needed migrating.
"""

import argparse
import shutil
import sys
from math import isclose
from pathlib import Path

from lxml import etree

# Messages of the files the run refused or could not process; decides the exit
# code.
errors: list[str] = []

# Relative to a <medium> element, so that a project file with several media is
# decided one medium at a time.
SOLID_COEFFICIENT_XPATH_IN_MEDIUM = (
    "phases/phase[type='Solid']/properties/"
    "property[name='thermal_osmosis_coefficient']"
)
MEDIUM_PROPERTY_XPATH_IN_MEDIUM = (
    "properties/property["
    "name='thermal_osmosis_coefficient' or name='thermal_osmosis_permeability']"
)


def find_prj_files(root_dir: Path) -> list[Path]:
    """Find all .prj files in subdirectories."""
    return list(root_dir.rglob("*.prj"))


def _fail(message: str) -> bool:
    """Print an ERROR message prefixed for the per-file log and return False.

    Counted, so that a run which refused a file exits non-zero even though
    refusing is, per file, a normal outcome.
    """
    errors.append(message)
    print(f"  ERROR: {message}")
    return False


def _components(element, what: str, prj_path: Path) -> list[float] | None:
    """
    Return the <value> of a Constant property as a list of float components.

    Returns None (after logging) if the property is not of type Constant or its
    value cannot be parsed.
    """
    property_type = element.findtext("type")
    if property_type != "Constant":
        _fail(
            f"Cannot convert: {what} has type '{property_type}', "
            f"only 'Constant' is supported, in {prj_path}"
        )
        return None

    value_str = element.findtext("value")
    if value_str is None:
        _fail(f"Cannot convert: {what} has no <value> in {prj_path}")
        return None

    try:
        return [float(v) for v in value_str.split()]
    except ValueError:
        _fail(
            f"Cannot convert: cannot parse {what} value "
            f"'{value_str.strip()}' in {prj_path}"
        )
        return None


def _diagonal(components: list[float], what: str, prj_path: Path) -> list[float] | None:
    """
    Return the diagonal of a scalar, 2x2, or 3x3 Constant value.

    Returns None (after logging) for an unexpected number of components or a
    non-diagonal tensor: the conversion inverts kT = epsilon_T k / mu per
    component, which reproduces the tensor relation only when both tensors are
    diagonal.
    """
    number_of_components = len(components)
    if number_of_components == 1:
        return components

    dimension = {4: 2, 9: 3}.get(number_of_components)
    if dimension is None:
        _fail(
            f"Cannot convert: {what} has {number_of_components} components, "
            f"expected 1, 4, or 9, in {prj_path}"
        )
        return None

    if any(
        components[row * dimension + column] != 0.0
        for row in range(dimension)
        for column in range(dimension)
        if row != column
    ):
        _fail(f"Cannot convert: {what} is not diagonal in {prj_path}")
        return None

    return [components[i * dimension + i] for i in range(dimension)]


def _scalar_epsilon_T(diagonal: list[float], prj_path: Path) -> str | None:
    """
    Format the per-component epsilon_T as the scalar OGS Constant <value>.

    thermal_osmosis_permeability is a scalar property, so a conversion whose
    components differ has no equivalent single value and is rejected instead of
    silently keeping one of them.
    """
    first = diagonal[0]
    if any(not isclose(component, first, rel_tol=1e-12) for component in diagonal[1:]):
        _fail(
            f"Cannot convert: thermal_osmosis_permeability is a scalar, but "
            f"kT * mu / k differs per component, epsilon_T={diagonal} in "
            f"{prj_path}"
        )
        return None
    return repr(first)


def _thermo_osmotic_permeability_value(
    kT_element, medium, prj_path: Path
) -> str | None:
    """
    Return the <value> of the epsilon_T property equivalent to the given
    thermal_osmosis_coefficient element.

    epsilon_T = kT * (mu / k) is evaluated per diagonal component, with the
    intrinsic permeability k and the AqueousLiquid viscosity mu taken from
    `medium`. Returns None (after logging) if the medium does not provide
    exactly one of each, or if a value cannot be converted.
    """
    permeability_elements = medium.xpath("properties/property[name='permeability']")
    if len(permeability_elements) != 1:
        _fail(
            f"Cannot convert: expected exactly one permeability property in "
            f"the medium, found {len(permeability_elements)} in {prj_path}"
        )
        return None

    viscosity_elements = medium.xpath(
        "phases/phase[type='AqueousLiquid']/properties/property[name='viscosity']"
    )
    if len(viscosity_elements) != 1:
        _fail(
            f"Cannot convert: expected exactly one AqueousLiquid.viscosity "
            f"property in the medium, found {len(viscosity_elements)} in {prj_path}"
        )
        return None

    kT_components = _components(kT_element, "thermal_osmosis_coefficient", prj_path)
    permeability_components = _components(
        permeability_elements[0], "permeability", prj_path
    )
    viscosity_components = _components(viscosity_elements[0], "viscosity", prj_path)
    if (
        kT_components is None
        or permeability_components is None
        or viscosity_components is None
    ):
        return None

    if len(viscosity_components) != 1:
        _fail(
            f"Cannot convert: viscosity must be a scalar, found "
            f"{len(viscosity_components)} components in {prj_path}"
        )
        return None
    viscosity = viscosity_components[0]

    kT_diagonal = _diagonal(kT_components, "thermal_osmosis_coefficient", prj_path)
    permeability_diagonal = _diagonal(permeability_components, "permeability", prj_path)
    if kT_diagonal is None or permeability_diagonal is None:
        return None

    # A scalar stands for an isotropic tensor, so broadcast it to the other
    # value's dimension; differing tensor dimensions are a project file error.
    if len(kT_diagonal) == 1:
        kT_diagonal = kT_diagonal * len(permeability_diagonal)
    if len(permeability_diagonal) == 1:
        permeability_diagonal = permeability_diagonal * len(kT_diagonal)
    if len(kT_diagonal) != len(permeability_diagonal):
        _fail(
            f"Cannot convert: thermal_osmosis_coefficient and permeability "
            f"have different dimensions in {prj_path}"
        )
        return None

    if viscosity <= 0 or any(k <= 0 for k in permeability_diagonal):
        _fail(
            f"Cannot convert: permeability and viscosity must be > 0, got "
            f"k={permeability_diagonal}, mu={viscosity} in {prj_path}"
        )
        return None

    epsilon_T_diagonal = [
        kT * (viscosity / k)
        for kT, k in zip(kT_diagonal, permeability_diagonal, strict=True)
    ]
    print(
        f"  Converting: kT={kT_diagonal}, k={permeability_diagonal}, "
        f"mu={viscosity} -> epsilon_T={epsilon_T_diagonal}"
    )
    return _scalar_epsilon_T(epsilon_T_diagonal, prj_path)


def move_thermal_osmosis_coefficient(pf, prj_path: Path, convert: bool = False):
    """
    Move the thermo-osmosis property from the solid phase to the medium level,
    optionally converting it from thermal_osmosis_coefficient to
    thermal_osmosis_permeability.

    Args:
        pf: root Element of the parsed project file
        prj_path: path of the project file `pf` was parsed from
        convert: If True, rename the property to thermal_osmosis_permeability
                 and apply the conversion epsilon_T = kT * (mu / k)

    Returns:
        bool: True if the tree was modified
    """
    media = pf.xpath("//medium")
    if not media:
        print("  - No changes needed")
        return False

    # Every medium is decided on its own: a project file may define one medium
    # that still carries the property on its solid phase next to another that
    # was migrated already, and a file-wide count would refuse both.
    modified = False
    for medium in media:
        if _move_in_medium(medium, prj_path, convert):
            modified = True

    if not modified:
        print("  - No changes needed")
    return modified


def _move_in_medium(medium, prj_path: Path, convert: bool) -> bool:
    """Move the thermo-osmosis property of a single <medium>, see
    move_thermal_osmosis_coefficient()."""
    solid_elements = medium.xpath(SOLID_COEFFICIENT_XPATH_IN_MEDIUM)
    medium_elements = medium.xpath(MEDIUM_PROPERTY_XPATH_IN_MEDIUM)

    if medium_elements and solid_elements:
        return _fail(
            f"Ambiguous: thermal_osmosis_coefficient is defined on the solid "
            f"phase and a thermo-osmosis property already exists on the medium "
            f"level in {prj_path}"
        )
    if medium_elements:
        # Already migrated. Re-running --convert here would apply
        # epsilon_T = kT * (mu / k) a second time and silently corrupt the
        # value, so skip instead.
        print("  - Already migrated, no changes needed")
        return False
    if not solid_elements:
        return False
    if len(solid_elements) > 1:
        return _fail(
            f"Expected at most one solid-phase thermal_osmosis_coefficient "
            f"property per medium, found {len(solid_elements)} in {prj_path}"
        )

    kT_element = solid_elements[0]

    # The lookups below (and the eventual insertion) stay inside this medium: an
    # unscoped search or insertion could silently touch another medium's
    # permeability, another phase's viscosity, or another medium's <properties>
    # block.
    def target_location():
        """The medium's <properties>, created if the medium has none.

        <properties> is optional under <medium>, so a medium keeping all its
        properties on the phases has none to move the property into. Called
        only once the move is going ahead, so that a refused medium is left
        exactly as it was found.
        """
        medium_properties = medium.xpath("properties")
        return (
            medium_properties[0]
            if medium_properties
            else etree.SubElement(medium, "properties")
        )

    if convert:
        epsilon_T_value = _thermo_osmotic_permeability_value(
            kT_element, medium, prj_path
        )
        if epsilon_T_value is None:
            return False
        target_location().append(
            etree.XML(
                "<property>"
                "<name>thermal_osmosis_permeability</name>"
                "<type>Constant</type>"
                f"<value>{epsilon_T_value}</value>"
                "</property>"
            )
        )
        kT_element.getparent().remove(kT_element)
        print(
            "  OK: Converted thermal_osmosis_coefficient to "
            "thermal_osmosis_permeability on the medium level"
        )
    else:
        # Without --convert the property's name and physical meaning do not
        # change, only its location. Re-parenting the element itself (lxml
        # detaches it from the solid phase) keeps the value verbatim, including
        # tensor values and non-Constant property types.
        target_location().append(kT_element)
        print("  OK: Moved thermal_osmosis_coefficient from solid to medium level")

    return True


def process_prj_file(
    prj_path: Path,
    dry_run: bool = False,
    backup: bool = False,
    convert: bool = False,
) -> bool:
    """
    Process a single PRJ file using lxml.etree.

    Args:
        prj_path: Path to PRJ file
        dry_run: If True, don't write changes
        backup: If True, create backup before modification
        convert: If True, apply conversion epsilon_T = kT * (mu / k)

    Returns:
        bool: True if file was modified
    """
    try:
        print(f"Processing: {prj_path}")

        pf = etree.parse(prj_path)
        pf_root = pf.getroot()

        status = move_thermal_osmosis_coefficient(pf_root, prj_path, convert)
        if status and not dry_run:
            # The backup is taken here, once the file is known to change, so
            # that a tree-wide run does not litter unmodified projects with
            # .prj.backup files.
            if backup:
                backup_path = prj_path.with_suffix(".prj.backup")
                shutil.copy2(prj_path, backup_path)
                print(f"  Created backup: {backup_path}")
            # Write in the encoding the file was read in. lxml defaults to
            # ASCII, which escapes every non-ASCII character as a numeric
            # character reference -- inside a comment, where such references
            # are not markup, that mangles the text permanently.
            etree.ElementTree(pf_root).write(
                prj_path,
                xml_declaration=True,
                encoding=pf.docinfo.encoding or "UTF-8",
            )
        return status

    except Exception as e:
        return _fail(f"{prj_path}: {e}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Move the thermo-osmosis property of every project file below "
            "ROOT_DIR from the solid phase to the medium level, optionally "
            "converting thermal_osmosis_coefficient to the equivalent scalar "
            "thermal_osmosis_permeability. Exits 1 if any file was refused or "
            "could not be parsed."
        )
    )
    parser.add_argument(
        "root_dir",
        nargs="?",
        default=Path.cwd(),
        type=Path,
        metavar="ROOT_DIR",
        help="directory searched recursively for .prj files (default: cwd)",
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="report changes without writing them"
    )
    parser.add_argument(
        "--backup",
        action="store_true",
        help="write a .prj.backup next to every file that changes",
    )
    parser.add_argument(
        "--convert",
        action="store_true",
        help=(
            "also convert the value: epsilon_T = kT * mu / k, using the "
            "medium's permeability and the aqueous liquid phase's viscosity"
        ),
    )
    args = parser.parse_args()
    root_dir = args.root_dir
    dry_run, backup, convert = args.dry_run, args.backup, args.convert

    print(f"Searching for PRJ files in: {root_dir}")
    print(f"Dry run: {dry_run}, Backup: {backup}, Convert: {convert}")
    print("-" * 60)

    # Find all PRJ files
    prj_files = find_prj_files(root_dir)
    print(f"Found {len(prj_files)} PRJ files")

    if not prj_files:
        print("No PRJ files found.")
        return

    # Process each file
    modified_count = 0
    for prj_path in prj_files:
        if process_prj_file(prj_path, dry_run, backup, convert):
            modified_count += 1

    print("\n" + "=" * 60)
    print(f"Summary: {modified_count} file(s) modified")
    if dry_run:
        print("(Dry run: no actual changes made)")
    if errors:
        print(f"{len(errors)} file(s) could not be migrated, see the ERROR lines")
        sys.exit(1)


if __name__ == "__main__":
    main()
