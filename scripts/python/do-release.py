import re
import subprocess
from pathlib import Path


def increment_ver(version):
    version = version.split(".")
    version[2] = str(int(version[2]) + 1)
    return ".".join(version)


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

# Modify changelog
new_version = ""
changelog_content = ""
with (source_path / "CHANGELOG.md").open() as f:
    changelog_content = f.read()

new_version_match = re.search(
    r"\[Please see the wiki-page\]\(.*Release-notes-([0-9]+\.[0-9]+\.[0-9]+)\)",
    changelog_content,
)
new_version = new_version_match[1]

print(f"New version: {new_version}")

changelog_new_content = re.sub(
    r"\[Please see the wiki-page\]\(.*Release-notes-[0-9]+\.[0-9]+\.[0-9]+\)\n\n----",
    f"""[Please see the wiki-page](https://gitlab.opengeosys.org/ogs/ogs/-/wikis/Release-notes-{increment_ver(new_version)})

## {new_version}

[Changelog for OpenGeoSys {new_version}](https://gitlab.opengeosys.org/ogs/ogs/-/wikis/Release-notes-{new_version})

----""",
    changelog_content,
)

with (source_path / "CHANGELOG.md").open("w") as f:
    f.write(changelog_new_content)

# Create web release page
subprocess.run(
    f"hugo new releases/{new_version}.md",
    shell=True,
    check=True,
    cwd=source_path / "web",
)

path = Path(source_path / "Documentation" / "mainpage.dox.in")
text = path.read_text()
text = text.replace(
    " * The documentation for OGS releases can be found here:\n *\n",
    f" * The documentation for OGS releases can be found here:\n *\n * - https://doxygen.opengeosys.org/{new_version}\n",
)
path.write_text(text)

path = Path(source_path / "README.md")
text = path.read_text()
text = text.replace(
    f"https://doxygen.opengeosys.org/{current_version}",
    f"https://doxygen.opengeosys.org/{new_version}",
)
path.write_text(text)

new_version_dash = new_version.replace(".", "-")
with (source_path / "scripts" / "doc" / "_redirects").open("a") as f:
    f.write(
        f"/{new_version}/* https://ogs-doxygen-{new_version_dash}.netlify.app/:splat 200!\n"
    )

print(
    f"""Run the following and update CITATION.cff:

git shortlog -sne HEAD...{current_version}

Then check diff and run:

git commit -m "{new_version}"
git tag -s -m "OpenGeoSys {new_version}" {new_version}
git push --tags"""
)
