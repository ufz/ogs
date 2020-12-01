# Inserts hashes of description + file path as GitLab code quality fingerprint
import hashlib
import json
import sys

data = None
with open(sys.argv[1]) as json_file:
    data = json.load(json_file)

for entry in data:
    desc = entry["description"]
    path = entry["location"]["path"]
    hash = hashlib.sha256((desc + path).encode("utf-8")).hexdigest()
    entry["fingerprint"] = hash

with open(sys.argv[1], "w") as outfile:
    json.dump(data, outfile)

print("Added cppcheck fingerprints to {}.".format(sys.argv[1]))
