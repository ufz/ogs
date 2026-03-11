# /// script
# dependencies = [
#   "lxml",
#   "numpy",
#   "pandas",
#   "typing_extensions",
# ]
# ///

import argparse
import itertools
from pathlib import Path

from FindFeatures import DIFF_FILE_FALLBACK_BASES, get_feature_dict, load_json_features
from utils import get_xml_files


def test_feature_matrix(path: Path, path_to_xml_files: Path) -> None:
    """Test the feature matrix."""
    xmls = get_xml_files(
        path_to_xml_files,
        process_diff_files=True,
        fallback_bases=DIFF_FILE_FALLBACK_BASES,
    )

    features = load_json_features(path)
    all_features = get_all_features(features)

    feature_dict = get_feature_dict(path, xml_files=xmls)

    missing_features = set(feature_dict.keys()) - set(all_features)

    # TODO: Remove the exception for the submesh residuum output once a test is implemented
    assert len(missing_features) == 0 or missing_features == {
        "Submesh Residuum Output: XDMF"
    }, (
        "The following features are present in the feature dictionary but could not be found in the project files: "
        + str(missing_features)
    )  # Checks whether all features appear in the feature dictionary

    assert len(features) == len(
        xmls
    ), "The following parseable xml files are missing in the feature matrix: " + str(
        [
            xml.docinfo.URL
            for xml in xmls
            if xml.docinfo.URL
            not in ["Tests/Data/" + feat["file"] for feat in features]
        ]
    )  # Checks whether all parsed xml files appear in the feature martrix

    assert len(all_features) >= 132  # Checks minimum number of features

    for feat in features:
        assert (
            len(feat["features"]) >= 3
        )  # Checks minimum number of features per xml file
        assert len(feat["features"]) == len(
            feat["lines"]
        )  # Check whether there is a line for each feature
        assert all(
            key in feat for key in ["file", "features", "lines", "feature_coverage"]
        )  # Check whether all keys are there

        # Check for feature coverage. If the feature coverage lies below that indicates problems with the feature detection and hence should raise an error.
        assert feat["feature_coverage"] > 0.94, (
            "The following file has a critically low feature coverage: " + feat["file"]
        )


def get_all_features(features: list[dict]) -> list[str]:
    """Returns a list of all features."""
    return list(
        set(itertools.chain.from_iterable([feat["features"] for feat in features]))
    )


if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Test generated feature matrix.")
    parser.add_argument(
        "path",
        type=Path,
        metavar="PATH",
        nargs="?",
        default="./build/DocAux/features.json",
        help="Path to the feature matrix.",
    )
    parser.add_argument(
        "--path_prj",
        type=Path,
        help="Path to the folder, where the project files are stored. Usually something like ./Tests/Data",
        default="./Tests/Data",
    )

    args = parser.parse_args()

    test_feature_matrix(args.path, args.path_prj)
