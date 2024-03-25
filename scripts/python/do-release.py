import re
import subprocess
from pathlib import Path


def increment_ver(version):
    version = version.split(".")
    version[2] = str(int(version[2]) + 1)
    return ".".join(version)


script_path = Path(__file__).resolve()
source_path = script_path.parent.parent.parent.resolve()

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
    f" * The documentation for OGS releases can be found here:\n *\n * - https://doxygen.opengeosys.org/v{new_version}\n",
)
path.write_text(text)
