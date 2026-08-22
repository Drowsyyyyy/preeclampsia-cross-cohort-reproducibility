"""Download the public source files used by the project.

Existing files are retained when their SHA-256 hash agrees with the recorded
manifest. A mismatching file is never silently accepted.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import requests
import yaml


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
URLS = ROOT / "config" / "source_urls.yaml"
MANIFEST = ROOT / "data" / "processed" / "data_manifest.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def expected_hashes() -> dict[str, str]:
    if not MANIFEST.exists():
        return {}
    payload = json.loads(MANIFEST.read_text(encoding="utf-8"))
    return {Path(item["path"]).name: item["sha256"] for item in payload["files"]}


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    urls = yaml.safe_load(URLS.read_text(encoding="utf-8"))
    hashes = expected_hashes()

    for name, url in urls.items():
        destination = RAW / name
        if destination.exists():
            observed = sha256(destination)
            expected = hashes.get(name)
            if expected and observed != expected:
                raise RuntimeError(
                    f"Existing file hash mismatch for {name}: {observed} != {expected}. "
                    "Remove or quarantine the file before retrying."
                )
            print(f"verified existing {name} ({destination.stat().st_size:,} bytes)")
            continue

        partial = destination.with_suffix(destination.suffix + ".part")
        print(f"downloading {name}")
        with requests.get(url, stream=True, timeout=120) as response:
            response.raise_for_status()
            with partial.open("wb") as handle:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        handle.write(chunk)
        partial.replace(destination)

        observed = sha256(destination)
        expected = hashes.get(name)
        if expected and observed != expected:
            raise RuntimeError(f"Downloaded file hash mismatch for {name}: {observed} != {expected}")
        print(f"downloaded {name} ({destination.stat().st_size:,} bytes; sha256={observed})")


if __name__ == "__main__":
    main()
