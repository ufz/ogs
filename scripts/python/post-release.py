import json
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

script_path = Path(__file__).resolve()
source_path = script_path.parent.parent.parent.resolve()
software_bib_path = source_path / "web/content/publications/software.bib"


@dataclass(frozen=True)
class SoftwareRelease:
    key_prefix: str
    concept_doi: str

    @property
    def concept_record_id(self) -> str:
        return self.concept_doi.rsplit(".", 1)[-1]

    @property
    def latest_key(self) -> str:
        return f"{self.key_prefix}:latest"


@dataclass(frozen=True)
class ZenodoRelease:
    record_id: str
    version: str
    bibtex: str


software_releases = (
    SoftwareRelease("ogs", "10.5281/zenodo.591265"),
    SoftwareRelease("ogstools", "10.5281/zenodo.10912873"),
)


def curl(*args: str) -> str:
    result = subprocess.run(
        ["curl", "--fail", "-L", "-s", *args],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout


def fetch_latest_zenodo_release(release: SoftwareRelease) -> ZenodoRelease:
    api_url = (
        f"https://zenodo.org/api/records/{release.concept_record_id}/versions/latest"
    )
    metadata = json.loads(curl(api_url))
    record_id = str(metadata["id"])
    version = metadata["metadata"]["version"]
    bibtex = curl(f"https://zenodo.org/records/{record_id}/export/bibtex")
    return ZenodoRelease(
        record_id=record_id,
        version=version,
        bibtex=rewrite_bibtex_key(bibtex.strip(), release.latest_key),
    )


def rewrite_bibtex_key(entry: str, key: str) -> str:
    return re.sub(r"^@([A-Za-z]+)\{[^,\n]+,", rf"@\1{{{key},", entry, count=1)


def bibtex_key(entry: str) -> str:
    match = re.match(r"^@[A-Za-z]+\{([^,\n]+),", entry.strip())
    if not match:
        msg = f"Could not determine BibTeX key for entry:\n{entry}"
        raise RuntimeError(msg)
    return match.group(1)


def bibtex_version(entry: str) -> str | None:
    match = re.search(r"\bversion\s*=\s*\{([^{}]+)\}", entry)
    if not match:
        return None
    return match.group(1).strip()


def split_bibtex_entries(content: str) -> list[str]:
    entries = []
    cursor = 0
    while True:
        start = content.find("@", cursor)
        if start == -1:
            break

        brace_start = content.find("{", start)
        if brace_start == -1:
            msg = f"Malformed BibTeX entry starting at byte {start}"
            raise RuntimeError(msg)

        depth = 0
        end = None
        for index in range(brace_start, len(content)):
            char = content[index]
            if char == "{":
                depth += 1
            elif char == "}":
                depth -= 1
                if depth == 0:
                    end = index + 1
                    break

        if end is None:
            msg = f"Unclosed BibTeX entry starting at byte {start}"
            raise RuntimeError(msg)

        entries.append(content[start:end].strip())
        cursor = end

    return entries


def update_release_entries(
    entries: list[str], release: SoftwareRelease, latest: ZenodoRelease
) -> list[str]:
    keys = [bibtex_key(entry) for entry in entries]
    latest_index = keys.index(release.latest_key) if release.latest_key in keys else -1

    if latest_index == -1:
        insert_index = next(
            (
                index
                for index, key in enumerate(keys)
                if key.startswith(f"{release.key_prefix}:")
            ),
            len(entries),
        )
        return [*entries[:insert_index], latest.bibtex, *entries[insert_index:]]

    updated_entries = entries.copy()
    former_latest = entries[latest_index]
    former_version = bibtex_version(former_latest)
    latest_version = latest.version

    updated_entries[latest_index] = latest.bibtex

    if former_version and former_version != latest_version:
        former_version_key = f"{release.key_prefix}:{former_version}"
        if former_version_key not in keys:
            updated_entries.insert(
                latest_index + 1,
                rewrite_bibtex_key(former_latest, former_version_key),
            )

    return updated_entries


def update_software_bib(latest_releases: dict[str, ZenodoRelease]) -> None:
    entries = split_bibtex_entries(software_bib_path.read_text())

    for release in software_releases:
        entries = update_release_entries(
            entries, release, latest_releases[release.key_prefix]
        )

    software_bib_path.write_text("\n\n".join(entries) + "\n")


def fetch_gitlab_release_changelog(version: str) -> str:
    release_json = curl(
        f"https://gitlab.opengeosys.org/api/v4/projects/120/releases/{version}"
    )
    return json.loads(release_json)["description"]


def main() -> None:
    latest_releases = {}
    for release in software_releases:
        latest_releases[release.key_prefix] = fetch_latest_zenodo_release(release)
        print(
            f"{release.key_prefix} latest Zenodo record: "
            f"{latest_releases[release.key_prefix].record_id}"
        )

    print("Issuing scan on Software Heritage ...")
    subprocess.run(
        [
            "curl",
            "-X",
            "POST",
            "https://archive.softwareheritage.org/api/1/origin/save/git/url/https://gitlab.opengeosys.org/ogs/ogs.git/",
        ],
        check=True,
    )

    print("\n\nUpdating software BibTeX entries on website ...")
    update_software_bib(latest_releases)

    ogs_latest = latest_releases["ogs"]
    changelog = fetch_gitlab_release_changelog(ogs_latest.version)

    print(f"""Add to CITATION.cff:

```
identifiers:
    - type: doi
      value: 10.5281/zenodo.{ogs_latest.record_id}
      description: Zenodo DOI for {ogs_latest.version}
```

Commit and push to master and create Discourse release post:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.{ogs_latest.record_id}.svg)](https://doi.org/10.5281/zenodo.{ogs_latest.record_id})

We are happy to announce the release of **OpenGeoSys {ogs_latest.version}**!

## Links

- [Release page on www.opengeosys.org](https://www.opengeosys.org/stable/releases/{ogs_latest.version}/)
- [GitLab Release](https://gitlab.opengeosys.org/ogs/ogs/-/releases/{ogs_latest.version})
- [Zenodo Release](https://zenodo.org/record/{ogs_latest.record_id})
- [PyPI Release](https://pypi.org/project/ogs/{ogs_latest.version}/)
- Also available on [Conda](https://anaconda.org/channels/conda-forge/packages/ogs/overview) and [Guix](https://hpc.guix.info/package/ogs-serial)

## Highlights

{changelog}
""")


if __name__ == "__main__":
    main()
