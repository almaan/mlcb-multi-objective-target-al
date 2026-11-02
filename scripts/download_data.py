# scripts/download_data.py
"""Fetch (where possible) and verify the raw inputs used in the manuscript.

  DepMap Public 26Q1  CRISPRGeneEffect.csv, Model.csv   -> data/raw/depmap/
  STRING v12.0        9606.protein.{info,links}.v12.0   -> data/string_dl/
  CORUM 5.3           humanComplexes.txt                -> data/corum_dl/

STRING is downloaded automatically. DepMap's download portal sits behind a
browser check, so the two DepMap files must be downloaded by hand from
https://depmap.org/portal/data_page/?tab=allData (release "DepMap Public 26Q1").
CORUM is tried through its official API; if that fails, download the human
complexes (txt) from https://mips.helmholtz-muenchen.de/corum/ by hand.

Every file is checked against the MD5 of the exact copy used for the manuscript.
CORUM serves only its *current* release, so a newer CORUM will (correctly) fail
the check -- the run still works, but results may differ slightly.

Run from the repository root:  python scripts/download_data.py
"""
from __future__ import annotations
import hashlib
import sys
import urllib.request
from pathlib import Path

FILES = [
    # (local path, md5, url or None for manual download)
    ("data/raw/depmap/CRISPRGeneEffect.csv", "e4f75f92348388459c91401d20a9724e", None),
    ("data/raw/depmap/Model.csv", "a15d75dffcc5219111ca39598948df9a", None),
    ("data/string_dl/9606.protein.info.v12.0.txt.gz", "1ca9efcefcc68b5cc3b037a335cb870b",
     "https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz"),
    ("data/string_dl/9606.protein.links.v12.0.txt.gz", "5b7c91d02926e723ad30e81ab9d009a4",
     "https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz"),
    ("data/corum_dl/humanComplexes.txt", "5cc8f6a98de30907f99f80adafce704a",
     "https://mips.helmholtz-muenchen.de/fastapi-corum/public/file/download_current_file"
     "?file_id=human&file_format=txt"),
]

MANUAL = {
    "data/raw/depmap/CRISPRGeneEffect.csv":
        "download 'CRISPRGeneEffect.csv' (DepMap Public 26Q1) from "
        "https://depmap.org/portal/data_page/?tab=allData",
    "data/raw/depmap/Model.csv":
        "download 'Model.csv' (DepMap Public 26Q1) from "
        "https://depmap.org/portal/data_page/?tab=allData",
    "data/corum_dl/humanComplexes.txt":
        "download the human complexes (txt) from https://mips.helmholtz-muenchen.de/corum/ "
        "and save as humanComplexes.txt",
}


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest: Path) -> bool:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    try:
        with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as f:
            while chunk := r.read(1 << 20):
                f.write(chunk)
    except Exception as e:  # noqa: BLE001 -- report and fall back to manual
        print(f"  download failed: {e}")
        tmp.unlink(missing_ok=True)
        return False
    tmp.rename(dest)
    return True


def main() -> int:
    missing, mismatched = [], []
    for rel, want, url in FILES:
        p = Path(rel)
        print(f"{rel}")
        if not p.exists() and url:
            print(f"  fetching {url}")
            fetch(url, p)
        if not p.exists():
            print("  MISSING")
            missing.append(rel)
            continue
        got = md5(p)
        if got == want:
            print("  ok (md5 matches manuscript input)")
        else:
            print(f"  MD5 MISMATCH: got {got}, manuscript used {want}")
            mismatched.append(rel)
    if missing:
        print("\nManual downloads still needed:")
        for rel in missing:
            print(f"  - {rel}: {MANUAL.get(rel, 'see header of this script')}")
    if mismatched:
        print("\nFiles differ from the manuscript inputs (different release?):")
        for rel in mismatched:
            print(f"  - {rel}")
    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main())
