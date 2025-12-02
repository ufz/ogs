#!/usr/bin/env python3
import os
import re

import requests

GITLAB_URL = os.environ["CI_SERVER_URL"]
PROJECT_ID = os.environ["CI_PROJECT_ID"]
TOKEN = os.environ["GITLAB_TOKEN"]  # CI_JOB_TOKEN or personal/project token
KEEP = int(os.environ.get("KEEP_DEV_VERSIONS", "5"))

headers = {"PRIVATE-TOKEN": TOKEN}

DEV_PATTERN = re.compile(r".*\.dev\d+$")  # matches versions like 1.0.0.dev12


def get_packages():
    url = f"{GITLAB_URL}/api/v4/projects/{PROJECT_ID}/packages?package_type=pypi&per_page=100"
    r = requests.get(url, headers=headers)
    r.raise_for_status()
    return r.json()


def delete_package(package_id):
    url = f"{GITLAB_URL}/api/v4/projects/{PROJECT_ID}/packages/{package_id}"
    r = requests.delete(url, headers=headers)
    r.raise_for_status()


def main():
    packages = get_packages()

    # Filter only .dev versions
    dev_packages = [p for p in packages if DEV_PATTERN.match(p["version"])]

    print(f"Found {len(dev_packages)} '.dev' packages")

    # Sort newest → oldest
    dev_packages = sorted(dev_packages, key=lambda p: p["created_at"], reverse=True)

    # Select which to delete
    to_delete = dev_packages[KEEP:]

    print(
        f"Keeping {KEEP} dev versions. Deleting {len(to_delete)} older dev packages..."
    )

    for pkg in to_delete:
        print(f"Deleting {pkg['name']} {pkg['version']} (id {pkg['id']})")
        delete_package(pkg["id"])

    print("Cleanup complete.")


if __name__ == "__main__":
    main()
