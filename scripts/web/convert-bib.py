#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "bibtexparser==1.4.4",
# ]
# ///

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import sys
import unicodedata
from pathlib import Path

import bibtexparser
from bibtexparser.bibdatabase import BibDatabase
from bibtexparser.bparser import BibTexParser
from bibtexparser.bwriter import BibTexWriter
from bibtexparser.latexenc import latex_to_unicode

IMAGE_SUFFIXES = {".png", ".jpg", ".jpeg", ".webp", ".gif", ".svg"}
REQUIRED_FIELDS = ("title", "author", "year")
EXTRA_FIELDS = (
    "booktitle",
    "volume",
    "number",
    "issue",
    "pages",
    "isbn",
    "address",
    "note",
    "institution",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="convert-bib.py",
        description=(
            "Convert BibTeX entries to Hugo publication bundles. If --input is "
            "omitted, all .bib files in --section are converted."
        ),
    )
    parser.add_argument("--input", "-i", type=Path)
    parser.add_argument(
        "--section", "-s", type=Path, default=Path("content/publications")
    )
    parser.add_argument("--output-section", "-o", type=Path)
    return parser.parse_args()


def decode_latex_value(value: str) -> str:
    decoded = latex_to_unicode(value)
    decoded = decoded.replace("\\textemdash", "—").replace("\\textendash", "–")
    decoded = decoded.replace("---", "—").replace("--", "–")
    return re.sub(r"\s+", " ", decoded).strip()


def slugify(value: str) -> str:
    normalised = unicodedata.normalize("NFKD", value.lower())
    ascii_only = "".join(ch for ch in normalised if not unicodedata.combining(ch))
    slug = re.sub(r"[^a-z0-9]+", "-", ascii_only).strip("-")
    return slug[:80]


def parse_authors(author_field: str) -> list[str]:
    authors = []
    for raw_author in re.split(r"\s+and\s+", author_field, flags=re.IGNORECASE):
        author = raw_author.strip()
        if not author:
            continue
        if "," in author:
            last, first = [part.strip() for part in author.split(",", 1)]
            author = " ".join(part for part in (first, last) if part)
        authors.append(author)
    return authors


def toml_string(value: object) -> str:
    return json.dumps(str(value), ensure_ascii=False)


def yaml_multiline(key: str, value: str) -> list[str]:
    return [f"{key}: |-", *[f"  {line}" for line in str(value).splitlines()]]


def create_display_title(entry: dict[str, str]) -> str:
    title = entry["title"]
    if entry["ENTRYTYPE"] != "software":
        return title
    version = entry.get("version", "").strip()
    return f"{title} {version}".strip()


def is_image_file(path: Path) -> bool:
    return path.suffix.lower() in IMAGE_SUFFIXES


def find_image_for_citation_key(
    citation_key: str, source_site_root: Path
) -> Path | None:
    image_base_dir = source_site_root / "content" / "img" / "publications"
    if not image_base_dir.is_dir():
        return None

    citation_key_lower = citation_key.lower()
    return next(
        (
            entry
            for entry in image_base_dir.iterdir()
            if entry.is_file()
            and is_image_file(entry)
            and entry.stem.lower() == citation_key_lower
        ),
        None,
    )


def copy_image_to_bundle(
    image_absolute_path: Path, bundle_path: Path
) -> tuple[str, Path]:
    bundle_image_path = bundle_path / image_absolute_path.name
    bundle_path.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(image_absolute_path, bundle_image_path)
    return image_absolute_path.name, bundle_image_path


def remove_existing_generated_files(bundle_path: Path) -> None:
    if not bundle_path.is_dir():
        return
    for entry in bundle_path.iterdir():
        if entry.is_file() and (entry.name == "index.md" or is_image_file(entry)):
            entry.unlink()


def find_nearest_hugo_root(start_dir: Path) -> Path | None:
    current_dir = start_dir.resolve()
    while True:
        if (current_dir / "hugo.toml").exists():
            return current_dir
        if current_dir.parent == current_dir:
            return None
        current_dir = current_dir.parent


def find_content_root(start_path: Path) -> Path | None:
    current_path = start_path.resolve()
    while True:
        if current_path.name == "content":
            return current_path
        if current_path.parent == current_path:
            return None
        current_path = current_path.parent


def resolve_output_section_path(
    section_path: Path, output_section_arg: Path | None
) -> Path:
    if output_section_arg is not None:
        return output_section_arg.resolve()
    return Path.cwd() / "content" / section_path.name


def format_bibtex_entry(entry: dict[str, str]) -> str:
    writer = BibTexWriter()
    writer.write_common_strings = False
    database = BibDatabase()
    database.entries = [entry]
    return bibtexparser.dumps(database, writer).strip()


def create_front_matter(entry: dict[str, str], image_path: str | None) -> str:
    title = create_display_title(entry)
    authors = parse_authors(entry["author"])
    year = int(entry["year"])
    journal = entry.get("journal")
    publisher = entry.get("publisher")
    venue = journal or publisher
    abstract = entry.get("abstract")
    summary = (
        f"{', '.join(authors)} ({year}). {venue}."
        if venue
        else f"{', '.join(authors)} ({year})."
    )

    lines = [
        "---",
        f"title: {toml_string(title)}",
        f"linkTitle: {toml_string(title)}",
        f"date: {toml_string(f'{year}-01-01')}",
        "weight: 1",
        "authors:",
        *[f"  - {toml_string(author)}" for author in authors],
        f"year: {year}",
        f"summary: {toml_string(summary)}",
        f"citationKey: {toml_string(entry['ID'])}",
        f"publicationType: {toml_string(entry.get('ENTRYTYPE', ''))}",
    ]

    if abstract:
        lines.append(f"abstract: {toml_string(abstract)}")
    if image_path:
        lines.append(f"image: {toml_string(image_path)}")
    if journal:
        lines.append(f"journal: {toml_string(journal)}")
    if publisher:
        lines.append(f"publisher: {toml_string(publisher)}")
    if entry.get("doi"):
        lines.append(f"doi: {toml_string(entry['doi'])}")
    if entry.get("url"):
        lines.append(f"sourceUrl: {toml_string(entry['url'])}")

    for field in EXTRA_FIELDS:
        if entry.get(field):
            lines.append(f"{field}: {toml_string(entry[field])}")

    lines.extend(yaml_multiline("bibtex", format_bibtex_entry(entry)))
    lines.extend(["---", ""])
    return "\n".join(lines)


def parse_bib_file(input_path: Path) -> list[dict[str, str]]:
    parser = BibTexParser(common_strings=True)
    parser.ignore_nonstandard_types = False
    parser.homogenize_fields = False

    with input_path.open(encoding="utf-8") as handle:
        database = bibtexparser.load(handle, parser=parser)

    entries = []
    for raw_entry in database.entries:
        entry = dict(raw_entry)
        for key, value in list(entry.items()):
            if key in {"ENTRYTYPE", "ID"} or not isinstance(value, str):
                continue
            entry[key] = decode_latex_value(value)
        entries.append(entry)
    return entries


def convert_bib_file(
    input_path: Path,
    output_section_path: Path,
    source_site_root: Path,
) -> int:
    entries = parse_bib_file(input_path)
    if not entries:
        msg = f"No BibTeX entries found in {input_path}"
        raise RuntimeError(msg)

    converted = 0
    for entry in entries:
        missing = [field for field in REQUIRED_FIELDS if not entry.get(field)]
        if missing:
            print(
                f"Skipping {entry['ID']}: missing required fields {', '.join(missing)}",
                file=sys.stderr,
            )
            continue

        slug_base = slugify(entry["ID"]) or slugify(entry["title"])
        bundle_path = output_section_path / slug_base
        image_src = find_image_for_citation_key(entry["ID"], source_site_root)
        image_file_name = None

        bundle_path.mkdir(parents=True, exist_ok=True)
        remove_existing_generated_files(bundle_path)
        if image_src is not None:
            image_file_name, _ = copy_image_to_bundle(image_src, bundle_path)
        page_content = create_front_matter(entry, image_file_name)
        (bundle_path / "index.md").write_text(page_content, encoding="utf-8")
        converted += 1
        image_suffix = " with image" if image_file_name is not None else ""
        print(
            f"Converted {entry['ID']} -> {os.path.relpath(bundle_path, Path.cwd())}{image_suffix}"
        )

    if converted == 0:
        msg = "No entries were converted. See warnings above."
        raise RuntimeError(msg)

    print(f"Done. Converted {converted} publication(s).")
    return converted


def find_bib_files(root_dir: Path) -> list[Path]:
    if not root_dir.exists():
        return []
    return sorted(path for path in root_dir.rglob("*.bib") if path.is_file())


def main() -> int:
    args = parse_args()
    section_path = args.section.resolve()
    output_section_path = resolve_output_section_path(section_path, args.output_section)
    source_site_root = (section_path / ".." / "..").resolve()
    bib_files = [args.input.resolve()] if args.input else find_bib_files(section_path)

    if not bib_files:
        msg = f"No .bib files found. Use --input or place .bib files in {section_path}"
        raise RuntimeError(msg)

    total = 0
    for bib_file in bib_files:
        total += convert_bib_file(
            bib_file,
            output_section_path,
            source_site_root,
        )

    print(
        f"Processed {len(bib_files)} .bib file(s). Converted {total} publication(s) in total."
    )
    return 0


def run() -> int:
    try:
        return main()
    except Exception as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(run())
