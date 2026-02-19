"""
File that contains utility functions for the creation of the documentation and Feature Matrix generation.
"""

import re
import warnings
from pathlib import Path
from typing import Literal

from lxml import etree

__all__ = [
    "add_includes",
    "apply_diff_tag",
    "apply_diff_to_base_xml",
    "get_diff_files",
    "get_diff_type",
    "get_xml_file_paths",
    "get_xml_files",
    "get_xpath",
    "process_diff_file",
    "replace_include_with_fragment",
    "replace_tag",
    "xml_parser",
]


def add_includes(
    xml: etree.ElementTree, root_dir: Path | None = None
) -> etree.ElementTree:
    """Add include references to the XML tree by replacing <include> elements with their content.

    This function processes all <include> elements in the XML tree, reads the referenced
    files, and replaces each include element with the content from the included file.

    Args:
        xml: The XML tree to process
        root_dir: Optional root directory. When provided, include paths that resolve
            outside this directory are rejected. Relative paths with '..' are allowed
            as long as they remain within root_dir.

    Returns:
        The XML tree with include elements replaced by their content
    """
    for include in xml.xpath("//include"):
        include_file = include.attrib["file"]
        # Reject absolute paths
        if include_file.startswith("/"):
            warnings.warn(
                f"Invalid include path '{include_file}': skipping include.",
                stacklevel=2,
            )
            continue
        include_path = Path(str(xml.docinfo.URL)).parent / include_file  # type: ignore[arg-type]
        if root_dir is not None:
            try:
                resolved = include_path.resolve()
                root_resolved = Path(root_dir).resolve()
                if not resolved.is_relative_to(root_resolved):
                    warnings.warn(
                        f"Invalid include path '{include_file}': resolves outside root directory, skipping include.",
                        stacklevel=2,
                    )
                    continue
            except (OSError, ValueError) as e:
                warnings.warn(
                    f"Could not resolve include path '{include_file}': {e}, skipping include.",
                    stacklevel=2,
                )
                continue
        replace_include_with_fragment(include, str(include_path))

    etree.indent(xml)  # Pretty-print after includes
    return xml


def apply_add_tag(
    diff_element: etree.Element,
    target_element: etree.Element,
    diff_type: Literal["text", "tag", "attrib"],
    attrib_name: str | None = None,
) -> None:
    """Add content to a target XML element based on the diff specification.

    This function applies "add" operations from a diff element to a target element.
    The type of addition depends on the diff_type parameter:

    - "text": Replaces the target element's text content with the diff element's text
    - "attrib": Adds or updates an attribute on the target element
    - "tag": Appends all non-comment child elements from the diff element to the target

    Args:
        diff_element: The XML element containing the addition specification and content
        target_element: The XML element to which content will be added
        diff_type: The type of addition operation ("text", "tag", or "attrib")
        attrib_name: The name of the attribute to add/update (required when diff_type is "attrib")

    Note:
        For "tag" type additions, comment nodes from the diff element are filtered out
        and not added to the target element.
    """
    match diff_type:
        case "text":
            change_text(target_element, diff_element.text)
        case "attrib":
            change_attrib(target_element, attrib_name, diff_element.text)
        case "tag":
            for child_to_add in list(diff_element):
                if not isinstance(child_to_add, etree._Comment):
                    target_element.append(
                        etree.fromstring(etree.tostring(child_to_add))
                    )


def apply_remove_tag(
    diff_element: etree.Element,  # noqa: ARG001 - kept for API consistency
    target_element: etree.Element,
    diff_type: Literal["text", "tag", "attrib"],
    attrib_name: str | None = None,
) -> None:
    """Remove content from a target XML element based on the diff specification.

    This function applies "remove" operations from a diff element to a target element.
    The type of removal depends on the diff_type parameter:

    - "text": Clears the target element's text content (sets it to empty string)
    - "attrib": Deletes a specific attribute from the target element
    - "tag": Removes the entire target element from its parent

    Args:
        diff_element: The XML element containing the removal specification
        target_element: The XML element from which content will be removed
        diff_type: The type of removal operation ("text", "tag", or "attrib")
        attrib_name: The name of the attribute to delete (required when diff_type is "attrib")

    Note:
        For "tag" type removals, the entire element is removed from its parent.
        This operation is irreversible and will delete the element and all its children.
    """
    match diff_type:
        case "text":
            target_element.text = ""
        case "attrib":
            del target_element.attrib[attrib_name]
        case "tag":
            target_element.getparent().remove(target_element)


def apply_replace_tag(
    diff_element: etree.Element,
    target_element: etree.Element,
    diff_type: Literal["text", "tag", "attrib"],
    attrib_name: str | None = None,
) -> None:
    """Replace content in a target XML element based on the diff specification.

    This function applies "replace" operations from a diff element to a target element.
    The type of replacement depends on the diff_type parameter:

    - "text": Replaces the target element's text content with the diff element's text
    - "attrib": Replaces the value of a specific attribute on the target element
    - "tag": Replaces the entire target element with the contents of the diff element

    Args:
        diff_element: The XML element containing the replacement specification and content
        target_element: The XML element to be replaced or modified
        diff_type: The type of replacement operation ("text", "tag", or "attrib")
        attrib_name: The name of the attribute to replace (required when diff_type is "attrib")

    Note:
        For "tag" type replacements, the entire target element is replaced with the
        children of the diff element. If the target element has no parent, a warning
        is issued and the operation is skipped. The function includes error handling
        to catch and warn about any exceptions during the replacement process.
    """
    match diff_type:
        case "text":
            change_text(target_element, diff_element.text)
        case "tag":
            parent = target_element.getparent()

            # Validate that target has a parent
            if parent is None:
                warnings.warn(
                    f"Target element has no parent, cannot replace: {etree.tostring(target_element, encoding='unicode')[:100]}",
                    stacklevel=2,
                )
                return
            try:
                # Replace the target element with the replacement
                replace_tag(target_element, diff_element)
            except (ValueError, etree.LxmlError) as e:
                warnings.warn(
                    f"Failed to replace element {etree.tostring(target_element, encoding='unicode')[:100]} "
                    f"with replacement {list(diff_element)}. Error: {e!s}",
                    stacklevel=2,
                )
        case "attrib":
            change_attrib(target_element, attrib_name, diff_element.text)


def change_attrib(
    elem: etree.Element, attribute_name: str, attribute_value: str
) -> None:
    """Set an attribute value on an XML element.

    Args:
        elem: The XML element to modify
        attribute_name: The name of the attribute to set
        attribute_value: The value to assign to the attribute
    """
    elem.attrib[attribute_name] = attribute_value


def change_text(elem: etree.Element, text: str) -> None:
    """Set the text content of an XML element.

    Args:
        elem: The XML element to modify
        text: The text content to assign to the element
    """
    elem.text = text


def get_attr_name(elem, xpath: str) -> str:
    """Extract the attribute name from an XPath expression or element type attribute.

    If the XPath ends with an attribute selector (starts with @), returns that attribute name.
    If the element has a 'type' attribute, extracts the attribute name from it.
    Otherwise returns None.

    Args:
        elem: The XML element containing the diff specification
        xpath: The XPath expression pointing to the target

    Returns:
        The attribute name if found, or None if not applicable

    Examples:
        >>> get_attr_name(elem, "path/to/@attribute")  # Returns "attribute"
        >>> get_attr_name(elem_with_type, "path/to")  # Returns "type_value" if type="@type"
    """
    last_segment = xpath.rsplit("/", maxsplit=1)[-1]
    if last_segment.startswith("@"):
        return last_segment[1:]
    if "type" in elem.attrib:
        return elem.attrib["type"][1:]
    return None


def get_diff_files(path: Path) -> list[Path]:
    """Extract all diff XML file paths from a directory.

    Scans the given directory recursively for .xml files that contain
    'OpenGeoSysProjectDiff' in their content, indicating they are diff files.

    Args:
        path: The directory path to search

    Returns:
        List of paths to diff XML files found
    """
    path_list_xml = list(path.rglob("*.xml"))
    path_list = []
    for filepath in path_list_xml:
        with filepath.open("rb") as included_xml_file:
            if not any(b"OpenGeoSysProjectDiff" in line for line in included_xml_file):
                continue
        path_list.append(filepath)
    return path_list


def get_diff_type(elem: etree.Element, xpath: str) -> Literal["text", "tag", "attrib"]:
    """Determine the type of diff operation based on XPath and element attributes.

    Classifies the diff operation as one of three types:
    - "text": when modifying text content
    - "attrib": when modifying attributes
    - "tag": when modifying element structure

    The classification considers both the XPath expression and element attributes:
    - If XPath ends with "text()", it's a text operation
    - If XPath targets an attribute (starts with @) or element has a "type" attribute, it's an attribute operation
    - Otherwise it's a tag operation

    Args:
        elem: The XML element containing the diff specification
        xpath: The XPath expression pointing to the target

    Returns:
        The diff type as a literal string ("text", "tag", or "attrib")
    """
    if xpath.endswith("text()"):
        return "text"
    if xpath.rsplit("/", maxsplit=1)[-1].startswith("@") or "type" in elem.attrib:
        return "attrib"
    return "tag"


def get_xml_file_paths(path: Path) -> list[Path]:
    """Extract all project file paths from a directory.

    Scans the given directory recursively for .prj files, which are OGS project files.

    Args:
        path: The directory path to search

    Returns:
        List of paths to project (.prj) files found
    """
    return list((path).rglob("*.prj"))


def get_xml_files(
    path: Path,
    process_diff_files: bool = False,
    fallback_bases: dict[str, str] | None = None,
) -> list[etree.ElementTree]:
    """Load and parse XML project files from a directory.

    Scans the given directory for .prj files and parses them into XML trees.
    Optionally processes diff files to apply modifications to base project files.

    Args:
        path: The directory path to search for project files
        process_diff_files: If True, also process .xml diff files that modify
                          base project files. Default is False.
        fallback_bases: Optional mapping from xml file path (relative to path)
            to base prj file (relative to the xml file's directory). Used for diff
            files that intentionally omit the base_file attribute.

    Returns:
        List of parsed XML trees for all project files (and processed diff files if enabled)

    Note:
        When process_diff_files is True, the function:
        1. Parses all .prj files
        2. Finds and processes .xml diff files
        3. Applies diff changes to base files in two passes (before and after includes)
        4. Returns combined list of all processed XML trees
    """
    files = get_xml_file_paths(path=path)

    if not files:
        msg = "No project files found"
        raise ValueError(msg)

    # Parse the XML Files
    if process_diff_files:
        xml_files = [add_includes(xml_parser(file), root_dir=path) for file in files]

        diff_files = get_diff_files(path=path)
        files.extend(diff_files)
        xml_files.extend(
            [
                process_diff_file(file, root_dir=path, fallback_bases=fallback_bases)
                for file in diff_files
            ]
        )

    else:
        xml_files = [xml_parser(file) for file in files]

    return [xml for xml in xml_files if xml is not None]


def get_xpath(elem: etree.Element) -> str:
    """Extract XPath expression from an element's attributes.

    Looks for either 'msel' or 'sel' attributes in the element and returns
    the corresponding XPath expression. 'msel' takes precedence if both exist.

    Args:
        elem: The XML element containing the XPath selector

    Returns:
        The XPath expression as a string
    """
    if "msel" in elem.attrib:
        xpath_expression = elem.attrib["msel"]
    else:
        xpath_expression = elem.attrib["sel"]
    return xpath_expression


def process_diff_file(
    path: Path,
    root_dir: Path | None = None,
    fallback_bases: dict[str, str] | None = None,
) -> etree.ElementTree | None:
    """Process a diff file and return the modified XML tree.

    This function reads a diff file, applies the specified changes to its
    base file, and returns the resulting XML tree. The processing happens
    in two passes:
    1. First pass: Apply changes before include resolution
    2. Add includes to the base XML
    3. Second pass: Apply changes after include resolution

    Args:
        path: Path to the diff XML file
        root_dir: Optional root directory. When provided, base file paths that
            resolve outside this directory are rejected. Relative paths with '..'
            are allowed as long as they remain within root_dir.
        fallback_bases: Optional mapping from xml file path (relative to root_dir)
            to base prj file (relative to the xml file's directory). Used for diff
            files that intentionally omit the base_file attribute because they can
            be applied to multiple base files.

    Returns:
        The modified XML tree with all diff changes applied, or None if
        the file is invalid or has no base file specified
    """
    xml_diff = xml_parser(path)
    if xml_diff is None:
        return None

    if "base_file" not in xml_diff.getroot().attrib:
        if fallback_bases is not None and root_dir is not None:
            rel_path = str(path.relative_to(root_dir))
            if rel_path not in fallback_bases:
                return None
            base_file = fallback_bases[rel_path]
        else:
            return None
    else:
        base_file = xml_diff.getroot().attrib["base_file"]

    # Reject absolute paths; relative paths with '..' are allowed within root_dir
    if base_file.startswith("/"):
        msg = f"Invalid base file path '{base_file}' in {xml_diff.docinfo.URL}: absolute paths are not allowed."
        warnings.warn(msg, stacklevel=2)
        return None

    base_dir = Path(str(xml_diff.docinfo.URL)).parent  # type: ignore[arg-type]
    base_xml_path = base_dir / base_file

    # Ensure the resolved path is within root_dir (or base_dir when root_dir is not given)
    try:
        base_xml_path = base_xml_path.resolve()
        boundary_dir = (
            Path(root_dir).resolve() if root_dir is not None else base_dir.resolve()
        )
        if not base_xml_path.is_relative_to(boundary_dir):
            msg = f"Base file path '{base_file}' resolves outside the {'root' if root_dir is not None else 'base'} directory."
            warnings.warn(msg, stacklevel=2)
            return None
    except (OSError, ValueError) as e:
        msg = f"Failed to resolve base file path '{base_file}': {e}"
        warnings.warn(msg, stacklevel=2)
        return None

    base_xml = xml_parser(base_xml_path)
    base_xml = apply_diff_to_base_xml(xml_diff, base_xml, False)
    base_xml = add_includes(base_xml, root_dir=root_dir)
    base_xml = apply_diff_to_base_xml(xml_diff, base_xml, True)
    base_xml.docinfo.URL = xml_diff.docinfo.URL
    return base_xml


def apply_diff_to_base_xml(
    xml: etree.ElementTree, base_xml: etree.ElementTree, after_include: bool
) -> etree.ElementTree:
    """Process diff XML elements and apply changes to the base XML tree.

    This function processes elements from a diff XML file and applies the corresponding
    modifications to a base XML tree. It handles three types of operations: replace,
    add, and remove, based on the element tags in the diff file.

    The function filters elements based on the 'after_includes' attribute. Elements
    with 'after_includes' set to True are processed only after includes have been
    resolved, while elements with 'after_includes' set to False are processed before
    includes are resolved. This two-pass approach ensures proper ordering when
    includes may affect which elements are available for modification.

    Args:
        xml: The diff XML tree containing the changes to be applied
        base_xml: The base XML tree to which changes will be applied
        after_include: Boolean flag indicating whether to process elements that
                      come after include resolution (True) or before (False)

    Returns:
        The modified base XML tree with diff changes applied

    Note:
        - Elements with 'after_includes' attribute equal to the 'after_include' parameter
          are processed
        - Comment elements are skipped
        - For each matching element, the XPath selector is used to locate target elements
        - If no targets are found for a selector, a LookupError is raised
        - The function supports text, attribute, and tag-level modifications
    """

    for elem in xml.getroot():
        after_includes_attribute = elem.attrib.get("after_includes", "false")

        # Convert string to boolean for proper comparison
        after_includes_bool = after_includes_attribute.lower() in ("true", "1")

        if after_includes_bool == after_include and not isinstance(
            elem, etree._Comment
        ):
            xpath_expression = get_xpath(elem)
            diff_type = get_diff_type(elem, xpath_expression)
            attrib_name = get_attr_name(elem, xpath_expression)

            if diff_type == "text" or xpath_expression.rsplit("/", maxsplit=1)[
                -1
            ].startswith("@"):
                xpath_expression = xpath_expression.rsplit("/", maxsplit=1)[0]

            targets = base_xml.xpath(xpath_expression)

            # Validate that we have targets to replace
            if not targets:
                selector = elem.attrib.get("sel", elem.attrib.get("msel", "N/A"))
                msg = f"No elements found for selector '{selector}'. Cannot perform replacement."
                raise LookupError(msg)

            for target_element in targets:
                apply_diff_tag(target_element, elem, diff_type, attrib_name)

    return base_xml


def apply_diff_tag(
    target_element: etree.Element,
    diff_element: etree.Element,
    diff_type: Literal["tag", "attrib", "text"],
    attrib_name: str | None = None,
) -> None:
    """Apply a diff operation (replace/add/remove) to a target element.

    This function dispatches the appropriate diff operation based on the
    diff_element's tag. It handles three types of modifications:
    - "replace": Replaces the target element or its content
    - "add": Adds new content to the target element
    - "remove": Removes the target element or its content

    Args:
        target_element: The XML element to modify
        diff_element: The diff specification element containing the operation
        diff_type: The type of diff operation ("tag", "attrib", or "text")
        attrib_name: The attribute name to modify (used for "attrib" type)

    Note:
        For "tag" type operations, the entire element structure is modified.
        For "attrib" type operations, only the specified attribute is changed.
        For "text" type operations, only the element's text content is modified.
    """
    match diff_element.tag:
        # Perform replacement for each target element

        case "replace":
            apply_replace_tag(diff_element, target_element, diff_type, attrib_name)

        case "add":
            apply_add_tag(diff_element, target_element, diff_type, attrib_name)

        case "remove":
            apply_remove_tag(diff_element, target_element, diff_type, attrib_name)

        case _:
            tag = diff_element.tag
            msg = f"Unknown tag: {tag}"
            raise ValueError(msg)


def replace_include_with_fragment(include_elem, include_path: str) -> None:
    """Replace an <include> element with the XML fragment contained in include_path.

    The include file may contain multiple top-level elements. The function
    wraps the content in a dummy root element for parsing, then inserts the
    parsed children at the include element's position in the parent tree.

    Args:
        include_elem: The <include> element to be replaced
        include_path: Path to the file containing the XML fragment to include
    """

    # Parse fragment by wrapping in a dummy root
    include_path_obj = Path(include_path)
    with include_path_obj.open("rb") as f:
        raw = f.read()

    raw = re.sub(rb"^\s*<\?xml\b[^?]*(?:\?(?!>)[^?]*)?\?>", b"", raw)

    diff_element = etree.fromstring(
        b"<fake_root>" + raw + b"</fake_root>",
        parser=etree.XMLParser(remove_blank_text=True),
    )
    replace_tag(include_elem, diff_element)


def replace_tag(element: etree.Element, replacement_element: etree.Element) -> None:
    """Replace an element with the contents of a replacement element.

    Inserts all children of the replacement element at the position of the
    original element, then removes the original element.

    Args:
        element: The XML element to be replaced
        replacement_element: The element whose children will replace the original

    Raises:
        ValueError: If the element has no parent (is a root element)
    """
    parent = element.getparent()
    if parent is None:
        msg = "include element has no parent"
        raise ValueError(msg)

    # Insert fragment children at the include's position
    insert_at = parent.index(element)
    for child in list(replacement_element):
        parent.insert(insert_at, etree.fromstring(etree.tostring(child)))
        insert_at += 1

    # Remove the original include tag
    parent.remove(element)


def xml_parser(filedir: Path) -> etree.ElementTree:
    """Parse an XML file into an ElementTree object.

    Wraps etree.parse with error handling to gracefully manage files that
    cannot be parsed due to various issues (malformed XML, encoding problems, etc.).

    Args:
        filedir: Path to the XML file to parse

    Returns:
        The parsed XML ElementTree, or None if parsing fails
    """
    try:
        return etree.parse(filedir)
    except Exception:
        warnings.warn(f"Not able to parse: {filedir}.", stacklevel=2)
