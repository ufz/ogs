import numpy as np
import pandas as pd
from lxml import etree


class feature_matrix:
    """A feature matrix object will be created. The evaluate feature matrix will set all the elements in the right place on a matrix.
    These are elements of the class "feature_matrix_element". All the attributes of these elements will also be attributes of the elements
    of this class, but in matrix shape."""

    def __init__(self, feature_dict, xml_files):
        self.feature_dict = feature_dict
        self.xml_files = {file.docinfo.URL: file for file in xml_files}
        self.feature_matrix_elements = feature_matrix.evaluateFeatureDict(
            self.xml_files, self.feature_dict
        )

        # Set all attributes of the feature_matrix_entry as attributes of this object and put them in matrix form.
        attribute_names = vars(self.feature_matrix_elements.iloc[0, 0])
        for attribute_name in list(attribute_names.keys()):
            setattr(
                self,
                attribute_name,
                self.feature_matrix_elements.map(
                    lambda x, attr=attribute_name: getattr(x, attr)
                ),
            )

        self.xml_length = {
            index: feature_matrix.getXMLEndLine(
                self.xml_files[index].xpath("../OpenGeoSysProject")[0]
            )
            for index in self.xml_files
        }
        (
            self.code_coverage,
            self.lines_without_features,
        ) = feature_matrix.getCodeCoverage(self.lines, self.xml_length)

    @staticmethod
    def getCodeCoverage(lines, xml_lengths) -> dict:
        """returns dictionary of percentage of code that is covered by the detected features."""
        lines_new = feature_matrix.getFeatureLines(lines)

        coverages = {
            index: sum(line.right - line.left + 1 for line in lines_new[index])
            / xml_lengths[index]
            for index, _row in lines.iterrows()
        }

        lines_no_feature = {}
        for file in lines_new:
            lines_no_feature.update(
                {
                    file: feature_matrix.getLinesNotCoveredByInterval(
                        lines_new[file], xml_lengths[file]
                    )
                }
            )
            # lines[mat.lines.index[0]].sort
        return coverages, lines_no_feature

    @staticmethod
    def getFeatureLines(lines) -> list[pd.Interval]:
        intervals_out = {}
        for index, row in lines.iterrows():
            lines_new = []
            [lines_new.append(j) for i in row for j in i]
            intervals = feature_matrix.mergeOverlappingIntervals(lines_new)
            intervals_out.update({index: intervals})
        return intervals_out

    @staticmethod
    def evaluateFeatureDict(
        xml_files: list[etree.ElementTree], featureDict: dict
    ) -> pd.DataFrame:
        """
        Created a matrix with rows = len(featureDict) and cols = len(xml_files).
        Each row represents a evaluation of the featureDict, while each column represents a xmlFile. Afterwards the values of the feature dict will be evaluated for each
        xmlFile.
        """

        # Initialize matrix with files as index and features as columns
        feature_matrix = pd.DataFrame(
            np.empty([len(xml_files), len(featureDict)], dtype=feature_matrix_entry),
            index=[xml_files[file].docinfo.URL for file in xml_files],
            columns=featureDict.keys(),
        )

        # Fill the matrix by iterating over each file and feature
        file_keys = list(xml_files.keys())
        for file_idx in range(len(xml_files)):
            features = [
                featureDict[key](xml_files[file_keys[file_idx]]) for key in featureDict
            ]
            feature_matrix.iloc[file_idx] = features

        return feature_matrix

    def getLinesNotCoveredByInterval(
        lines: list[pd.Interval], endpoint: int
    ) -> list[pd.Interval]:
        not_covered = []
        if len(lines) == 1:
            if lines[0].left != 1:
                not_covered.append(pd.Interval(1, lines[0].left - 1))
            if lines[0].right != endpoint:
                not_covered.append(pd.Interval(lines[0].right + 1, endpoint))

        for l in range(len(lines) - 1):
            if lines[l].left != 1 and l == 0:
                not_covered.append(pd.Interval(1, lines[l].left - 1, closed="both"))
            not_covered.append(
                pd.Interval(lines[l].right + 1, lines[l + 1].left - 1, closed="both")
            )
            if l == len(lines) - 2 and lines[l + 1].right < endpoint:
                not_covered.append(
                    pd.Interval(lines[l + 1].right + 1, endpoint, closed="both")
                )
        return not_covered

    @staticmethod
    def getXMLEndLine(element: etree.ElementTree) -> int:
        """Gets the sourceline where the tag is closed, by converting to string and checking the number of lines."""
        return element.sourceline + (
            len(etree.tostring(element).strip().split(b"\n")) - 1
        )

    @staticmethod
    def mergeOverlappingIntervals(
        intervals: list[pd.Interval],
    ) -> list[pd.Interval]:
        """Will merge intervals from a list. All intervals that overlap or are 1 step apart from each other will be merged into a single larger interval. Will return a cleaned list of intervals."""
        k = 0
        while k < (len(intervals) - 1):
            l = k + 1
            while l < len(intervals):
                if pd.Interval.overlaps(intervals[k], intervals[l]):
                    intervals[k] = pd.Interval(
                        min(intervals[k].left, intervals[l].left),
                        max(intervals[k].right, intervals[l].right),
                        closed="both",
                    )
                    del intervals[l]
                else:
                    l += 1
            k += 1
        m = 0
        intervals.sort()
        while m < (len(intervals) - 1):
            if intervals[m + 1].left - intervals[m].right <= 1:
                intervals[m] = pd.Interval(
                    intervals[m].left, intervals[m + 1].right, closed="both"
                )
                del intervals[m + 1]
            else:
                m += 1

        return intervals


class feature_matrix_entry:
    """
    A class for the entries of a feature matrix. All the attributes that this class contains will be attributes of the feature matrix.
    One or several etree.ElemenTree objects can be put in one feature_class object. If none is given, this class will interpret it as
    if the given feature is not present in the respective file.
    """

    def __init__(self, elements=list[etree.ElementTree], line_type="range", lines=None):
        if len(elements) > 0:
            if lines is not None:
                self.lines = lines
            else:
                match line_type:
                    case "range":
                        # Add the complete range of the Element to the lines
                        self.lines = [
                            pd.Interval(
                                el.sourceline,
                                feature_matrix.getXMLEndLine(el),
                                closed="both",
                            )
                            for el in elements
                        ]
                    case "open and close":
                        # Only add the lines, where the tag is opened (eg. <parameter>) and closed (eg.</parameter>)
                        self.lines = []
                        [
                            self.lines.extend(
                                [
                                    pd.Interval(
                                        el.sourceline,
                                        el.sourceline,
                                        closed="both",
                                    ),
                                    pd.Interval(
                                        feature_matrix.getXMLEndLine(el),
                                        feature_matrix.getXMLEndLine(el),
                                        closed="both",
                                    ),
                                ]
                            )
                            for el in elements
                        ]
            self.has_feature = True
        else:
            self.lines = []
            self.has_feature = False
