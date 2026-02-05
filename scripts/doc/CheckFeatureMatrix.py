import argparse
import itertools
from pathlib import Path

from FindFeatures import load_json_features
from utils import get_xml_files
from FindFeatures import get_feature_dict


def test_feature_matrix(path: Path, path_to_xml_files: Path) -> None:
    """Test the feature matrix."""
    xmls = get_xml_files(path_to_xml_files, False)

    features = load_json_features(path)
    all_features = get_all_features(features)

    feature_dict = get_feature_dict(path, xml_files=xmls)

    assert len(feature_dict) == len(
        all_features
    ), "The following features are missing in the feature dictionary: " + str(
        set(feature_dict.keys()) - set(all_features)
    )  # Checks whether all features appear in the feature dictionary

    assert len(features) == len(
        xmls
    ), "There seem to be parseable xml files missing in the feature matrix."  # Checks whether all parsed xml files appear in the feature martrix

    assert len(all_features) >= 132  # Checks minimum number of features

    for feat in features:
        assert (
            len(feat["features"]) >= 3
        )  # Checks minimum number of features per xml file
        assert len(feat["features"]) == len(
            feat["lines"]
        )  # Check whether there is a line for each feature
        assert not any(
            [
                mandatory_key not in list(feat.keys())
                for mandatory_key in [
                    "file",
                    "features",
                    "lines",
                    "feature_coverage",
                ]
            ]
        )  # Check whether all keys are there


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
