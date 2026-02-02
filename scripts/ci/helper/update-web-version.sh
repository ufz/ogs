#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <version> <subdomain> <json-file>" >&2
  exit 1
fi

VERSION="$1"
SUBDOMAIN="$2"
FILE="$3"

tmp="$(mktemp)"

jq --arg v "$VERSION" --arg sub "$SUBDOMAIN" '
  map(del(.preferred))
  | [
      (.[] | select(.version == "master")),
      {
        version: $v,
        url: ("https://" + $sub + ".opengeosys.org/" + $v + "/"),
        preferred: true
      },
      (.[] | select(.version != "master"))
    ]
' "$FILE" > "$tmp"

mv "$tmp" "$FILE"
