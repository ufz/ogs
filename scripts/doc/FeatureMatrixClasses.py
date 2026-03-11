import uuid
import numpy as np
import pandas as pd
from typing import Literal
from lxml import etree


class FeatureMatrix:
    """A feature matrix object will be created. The evaluate feature matrix will set all the elements in the right place on a matrix.
    These are elements of the class "FeatureMatrixElement". All the attributes of these elements will also be attributes of the elements
    of this class, but in matrix shape."""

    def __init__(self, feature_dict: dict, xml_files: list[etree.ElementTree]):
        self.feature_dict = feature_dict
        self.xml_files = {file.docinfo.URL: file for file in xml_files}
        self.feature_matrix_elements = FeatureMatrix.evaluate_feature_dict(
            self.xml_files, self.feature_dict
        )

        # Set all attributes of the feature_matrix_entry as attributes of this object and put them in matrix form.
        attribute_names = vars(self.feature_matrix_elements.iloc[0, 0])
        for attribute_name in attribute_names:
            setattr(
                self,
                attribute_name,
                self.feature_matrix_elements.map(
                    lambda x, attr=attribute_name: getattr(x, attr)
                ),
            )

        self.xml_length = {
            index: FeatureMatrix.get_xml_endline(
                self.xml_files[index].xpath("../OpenGeoSysProject")[0]
            )
            for index in self.xml_files
        }
        (
            self.code_coverage,
            self.lines_without_features,
        ) = FeatureMatrix.get_code_coverage(self.lines, self.xml_length)

    @staticmethod
    def get_code_coverage(
        lines: list[pd.Interval], xml_lengths: dict[str, int]
    ) -> tuple[dict, dict]:
        """returns dictionary of percentage of code that is covered by the detected features."""
        lines_new = FeatureMatrix.get_feature_lines(lines)

        coverages = {
            index: sum(line.right - line.left + 1 for line in lines_new[index])
            / xml_lengths[index]
            for index, _row in lines.iterrows()
        }

        lines_no_feature = {}
        for file in lines_new:
            lines_no_feature.update(
                {
                    file: FeatureMatrix.get_lines_not_covered_by_intervals(
                        lines_new[file], xml_lengths[file]
                    )
                }
            )
            # lines[mat.lines.index[0]].sort
        return coverages, lines_no_feature

    @staticmethod
    def get_feature_lines(lines: list[pd.Interval]) -> list[pd.Interval]:
        intervals_out = {}
        for index, row in lines.iterrows():
            lines_new = []
            [lines_new.append(j) for i in row for j in i]
            intervals = FeatureMatrix.merge_overlapping_intervals(lines_new)
            intervals_out.update({index: intervals})
        return intervals_out

    @staticmethod
    def evaluate_feature_dict(
        xml_files: dict[str, etree.ElementTree], featureDict: dict
    ) -> pd.DataFrame:
        """
        Created a matrix with rows = len(featureDict) and cols = len(xml_files).
        Each row represents a evaluation of the featureDict, while each column represents a xmlFile. Afterwards the values of the feature dict will be evaluated for each
        xmlFile.
        """

        # Initialize matrix with files as index and features as columns
        feature_matrix = pd.DataFrame(
            np.empty([len(xml_files), len(featureDict)], dtype=FeatureMatrixEntry),
            index=[xml_files[file].docinfo.URL for file in xml_files],
            columns=featureDict.keys(),
        )

        # Fill the matrix by iterating over each file and feature
        for file_idx, (file_key, xml_file) in enumerate(xml_files.items()):
            features = [featureDict[key](xml_file) for key in featureDict]
            feature_matrix.iloc[file_idx] = features

        return feature_matrix

    @staticmethod
    def get_lines_not_covered_by_intervals(
        lines: list[pd.Interval], endpoint: int
    ) -> list[pd.Interval]:
        not_covered = []
        if len(lines) == 1:
            if lines[0].left != 1:
                not_covered.append(pd.Interval(1, lines[0].left - 1, closed="both"))
            if lines[0].right != endpoint:
                not_covered.append(
                    pd.Interval(lines[0].right + 1, endpoint, closed="both")
                )
            return not_covered
        elif len(lines) == 0:
            return [pd.Interval(1, endpoint, closed="both")]
        for line_index in range(len(lines) - 1):
            if lines[line_index].left != 1 and line_index == 0:
                not_covered.append(
                    pd.Interval(1, lines[line_index].left - 1, closed="both")
                )
            not_covered.append(
                pd.Interval(
                    lines[line_index].right + 1,
                    lines[line_index + 1].left - 1,
                    closed="both",
                )
            )
            if line_index == len(lines) - 2 and lines[line_index + 1].right < endpoint:
                not_covered.append(
                    pd.Interval(
                        lines[line_index + 1].right + 1, endpoint, closed="both"
                    )
                )
        return not_covered

    # Class-level offset cache to track XML declaration offsets per file
    _xml_offset_cache: dict[str, int] = {}

    @staticmethod
    def _get_xml_offset(root: etree.ElementTree) -> int:
        """Calculate the line offset due to XML declaration and trailing newline being stripped by lxml.

        lxml's etree.parse() strips:
        1. The XML declaration (<?xml ...?>)
        2. Any trailing newline

        We need to add an offset to match line numbers with the original file.
        """
        # Use file URL as cache key (not id(root)) because getroottree() returns
        # a new _ElementTree wrapper on each call, so id() values are ephemeral
        # and can be reused after garbage collection, poisoning the cache.
        cache_key = root.docinfo.URL or ""
        if cache_key in FeatureMatrix._xml_offset_cache:
            return FeatureMatrix._xml_offset_cache[cache_key]

        # Check if the original file had an XML declaration
        # docinfo.xml_version is not None if the original file had <?xml ...?>
        offset = 1 if root.docinfo.xml_version is not None else 0

        # Check if the original file had a trailing newline
        # by comparing the serialized length to the original file length
        original_file = root.docinfo.URL
        if original_file:
            try:
                with open(original_file, "rb") as f:
                    original_content = f.read()
                serialized_content = etree.tostring(root)
                # If original ends with newline but serialized doesn't, add 1
                if original_content.endswith(b"\n") and not serialized_content.endswith(
                    b"\n"
                ):
                    offset += 1
            except (OSError, IOError):
                pass  # If we can't read the file, just use the XML declaration offset

        FeatureMatrix._xml_offset_cache[cache_key] = offset
        return offset

    @staticmethod
    def _get_element_line(element: etree.ElementTree, is_closing: bool) -> int:
        """Find the line number of an element's opening or closing tag in the original file.

        Uses a UUID marker to uniquely identify the element since multiple elements
        can have the same tag name. The tree should already be formatted with
        etree.indent() before this is called.
        """
        # Get the root tree (already formatted by utils.add_includes)
        root = element.getroottree()

        # Calculate offset for XML declaration
        xml_offset = FeatureMatrix._get_xml_offset(root)

        # Create a unique marker to identify this specific element
        marker_id = str(uuid.uuid4())
        element.set("data-uuid", marker_id)

        try:
            # Serialize the tree
            s = etree.tostring(root, encoding="unicode")
            lines = s.split("\n")

            if is_closing:
                # For closing tag, find </tag> after the opening tag with our UUID
                found_opening = False
                for i, line in enumerate(lines):
                    if f'data-uuid="{marker_id}"' in line:
                        found_opening = True
                    if found_opening and f"</{element.tag}>" in line:
                        return i + 1 + xml_offset  # Convert to 1-indexed and add offset
                # Self-closing element - closing is same as opening
                if found_opening:
                    for i, line in enumerate(lines):
                        if f'data-uuid="{marker_id}"' in line:
                            return i + 1 + xml_offset
                return -1
            else:
                # For opening tag, find the line with our unique marker
                for i, line in enumerate(lines):
                    if f'data-uuid="{marker_id}"' in line:
                        return i + 1 + xml_offset  # Convert to 1-indexed and add offset
                return -1
        finally:
            # Remove the marker from the element
            if "data-uuid" in element.attrib:
                del element.attrib["data-uuid"]

    @staticmethod
    def get_xml_endline(element: etree.ElementTree) -> int:
        """Gets the line where the closing tag is located."""
        return FeatureMatrix._get_element_line(element, is_closing=True)

    @staticmethod
    def get_element_opening_line(element):
        """Gets the line where the opening tag is located."""
        return FeatureMatrix._get_element_line(element, is_closing=False)

    @staticmethod
    def merge_overlapping_intervals(
        intervals: list[pd.Interval],
    ) -> list[pd.Interval]:
        """Will merge intervals from a list. All intervals that overlap or are 1 step apart from each other will be merged into a single larger interval. Will return a cleaned list of intervals."""
        interval_index_left = 0
        while interval_index_left < (len(intervals) - 1):
            interval_index_right = interval_index_left + 1
            while interval_index_right < len(intervals):
                if pd.Interval.overlaps(
                    intervals[interval_index_left], intervals[interval_index_right]
                ):
                    intervals[interval_index_left] = pd.Interval(
                        min(
                            intervals[interval_index_left].left,
                            intervals[interval_index_right].left,
                        ),
                        max(
                            intervals[interval_index_left].right,
                            intervals[interval_index_right].right,
                        ),
                        closed="both",
                    )
                    del intervals[interval_index_right]
                else:
                    interval_index_right += 1
            interval_index_left += 1
        m = 0
        intervals.sort()
        while m < (len(intervals) - 1):
            if intervals[m + 1].left - intervals[m].right <= 1:
                intervals[m] = pd.Interval(
                    intervals[m].left,
                    max(intervals[m].right, intervals[m + 1].right),
                    closed="both",
                )
                del intervals[m + 1]
            else:
                m += 1

        return intervals


class FeatureMatrixEntry:
    """
    A class for the entries of a feature matrix. All the attributes that this class contains will be attributes of the feature matrix.
    One or several etree.ElemenTree objects can be put in one feature_class object. If none is given, this class will interpret it as
    if the given feature is not present in the respective file.
    """

    @staticmethod
    def _create_range_intervals(
        elements: list[etree.ElementTree],
    ) -> list[pd.Interval]:
        """Create intervals spanning the complete range of each element from opening to closing tag."""
        intervals = []
        for el in elements:
            open_line = FeatureMatrix.get_element_opening_line(el)
            close_line = FeatureMatrix.get_xml_endline(el)
            # Ensure left <= right to avoid ValueError in pd.Interval
            if open_line <= close_line:
                intervals.append(pd.Interval(open_line, close_line, closed="both"))
            else:
                # If opening line > closing line, use a single-point interval at the opening line
                # This handles edge cases where line detection fails
                intervals.append(pd.Interval(open_line, open_line, closed="both"))
        return intervals

    @staticmethod
    def _create_open_close_intervals(
        elements: list[etree.ElementTree],
    ) -> list[pd.Interval]:
        """Create two single-line intervals for each element: one at opening tag line and one at closing tag line."""
        intervals = []
        for el in elements:
            open_line = FeatureMatrix.get_element_opening_line(el)
            close_line = FeatureMatrix.get_xml_endline(el)
            intervals.extend(
                [
                    pd.Interval(open_line, open_line, closed="both"),
                    pd.Interval(close_line, close_line, closed="both"),
                ]
            )
        return intervals

    def __init__(
        self,
        elements: list[etree.ElementTree],
        line_type: Literal["range", "open and close"] = "range",
        lines: list[pd.Interval] | None = None,
    ):
        if len(elements) > 0:
            if lines is not None:
                self.lines = lines
            else:
                match line_type:
                    case "range":
                        self.lines = FeatureMatrixEntry._create_range_intervals(
                            elements
                        )
                    case "open and close":
                        self.lines = FeatureMatrixEntry._create_open_close_intervals(
                            elements
                        )
            self.has_feature = True
        else:
            self.lines = []
            self.has_feature = False
