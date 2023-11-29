# Inserts hashes of description + file path as GitLab code quality fingerprint
import hashlib
import json
import sys
from pathlib import Path

data = None
with Path(sys.argv[1]).open() as json_file:
    data = json.load(json_file)

for entry in data:
    desc = entry["description"]
    path = entry["location"]["path"]
    hash = hashlib.sha256((desc + path).encode("utf-8")).hexdigest()
    entry["fingerprint"] = hash

with Path(sys.argv[1]).open("w") as outfile:
    json.dump(data, outfile)

print(f"Added cppcheck fingerprints to {sys.argv[1]}.")
