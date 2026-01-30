import argparse
import subprocess
from pathlib import Path


def increment_ver(version):
    version = version.split(".")
    version[2] = str(int(version[2]) + 1)
    return ".".join(version)


def main():
    parser = argparse.ArgumentParser(description="Update release files for OpenGeoSys")
    parser.add_argument(
        "--new_version", required=True, help="The new version number (e.g., 6.5.6)"
    )
    args = parser.parse_args()

    script_path = Path(__file__).resolve()
    source_path = script_path.parent.parent.parent.resolve()

    git_describe = subprocess.run(
        "git describe --tags --abbrev=0",
        shell=True,
        check=True,
        text=True,
        capture_output=True,
    )
    current_version = git_describe.stdout.splitlines()[0]
    print(f"Current version: {current_version}, new version: {args.new_version}")

    # Modify changelog
    changelog_content = ""
    with (source_path / "CHANGELOG.md").open() as f:
        changelog_content = f.read()

    # Insert new version section after the dashes (----) in the changelog
    new_version_section = f"""
## {args.new_version}

[Changelog for OpenGeoSys {args.new_version}](https://gitlab.opengeosys.org/ogs/ogs/-/wikis/Release-notes-{args.new_version})"""

    # split into lines (keeping line endings)
    lines = changelog_content.splitlines(keepends=True)

    # insert after line 7 (1-based), i.e. index 7 in 0-based list
    insert_at = 7
    lines.insert(insert_at, new_version_section + "\n")

    changelog_new_content = "".join(lines)

    with (source_path / "CHANGELOG.md").open("w") as f:
        f.write(changelog_new_content)

    # Create web release page
    subprocess.run(
        f"hugo new releases/{args.new_version}.md",
        shell=True,
        check=True,
        cwd=source_path / "web",
    )

    path = Path(source_path / "Documentation" / "mainpage.dox.in")
    text = path.read_text()
    text = text.replace(
        " * The documentation for OGS releases can be found here:\n *\n",
        f" * The documentation for OGS releases can be found here:\n *\n * - https://doxygen.opengeosys.org/{args.new_version}\n",
    )
    path.write_text(text)

    path = Path(source_path / "README.md")
    text = path.read_text()
    text = text.replace(
        f"https://doxygen.opengeosys.org/{current_version}",
        f"https://doxygen.opengeosys.org/{args.new_version}",
    )
    path.write_text(text)

    print(
        f"""Run the following and update CITATION.cff:

git shortlog -sne HEAD...{current_version}

Then check diff and run:

git commit -m "{args.new_version}"
git tag -s -m "OpenGeoSys {args.new_version}" {args.new_version}
git push origin --tags
git push origin master"""
    )


if __name__ == "__main__":
    main()
