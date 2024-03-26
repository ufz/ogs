import re
import subprocess
from pathlib import Path

script_path = Path(__file__).resolve()
source_path = script_path.parent.parent.parent.resolve()

r = subprocess.run(
    "curl https://zenodo.org/doi/10.5281/zenodo.591265",
    shell=True,
    check=True,
    capture_output=True,
    text=True,
)

zenodo_record = re.search(r'"https://zenodo.org/records/([0-9]+)"', r.stdout).group(1)
print(f"Zenodo record: {zenodo_record}")

print("Issuing scan on Software Heritage ...")
r = subprocess.run(
    "curl -X POST https://archive.softwareheritage.org/api/1/origin/save/git/url/https://gitlab.opengeosys.org/ogs/ogs.git/",
    shell=True,
    check=True,
)

print("\n\nUpdating bibtex entry on website ...")
r = subprocess.run(
    f"curl https://zenodo.org/records/{zenodo_record}/export/bibtex",
    shell=True,
    check=True,
    capture_output=True,
    text=True,
)

bibtex = r.stdout

m = re.search(r"\{([0-9]+\.[0-9]+\.[0-9]+)\}", bibtex)
version = m.group(1)

bibtex = re.sub(r"@software\{.*,\n", f"@software{{ogs:{version},\n", bibtex)

publications_index_page = ""
with (source_path / "web/content/publications/_index.md").open() as f:
    publications_index_page = f.read()

publications_index_page = re.sub(
    r"(?s)```bibtex.*```", f"```bibtex\n{bibtex}\n```", publications_index_page
)

with (source_path / "web/content/publications/_index.md").open("w") as f:
    f.write(publications_index_page)

r = subprocess.run(
    f"curl -s https://gitlab.opengeosys.org/api/v4/projects/120/releases/{version} | jq '.description'",
    shell=True,
    check=True,
    text=True,
    capture_output=True,
)
changelog = r.stdout[1:-2].replace(r"\n", "\n").replace(r"\"", '"')

print(
    f"""Add to CITATION.cff:

```
identifiers:
    - type: doi
      value: 10.5281/zenodo.{zenodo_record}
      description: Zenodo DOI for {version}
```

Commit and push to master and create Zenodo release post:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.{zenodo_record}.svg)](https://doi.org/10.5281/zenodo.{zenodo_record})

We are happy to announce the release of **OpenGeoSys {version}**!

## Links

- [Release page on www.opengeosys.org](https://www.opengeosys.org/releases/{version}/)
- [GitLab Release](https://gitlab.opengeosys.org/ogs/ogs/-/releases/{version})
- [Zenodo Release](https://zenodo.org/record/{zenodo_record})
- [PyPI Release](https://pypi.org/project/ogs/{version}/)

## Highlights

{changelog}
"""
)
