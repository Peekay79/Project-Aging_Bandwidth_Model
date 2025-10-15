#!/usr/bin/env python3
"""
Download and cache real datasets for the mouse aging headroom analysis.

This script is safe to re-run; existing files are not re-downloaded unless --force is used.

Datasets:
1) Mouse Aging Proteome Atlas (TMT proteomics across multiple tissues and ages)
   - Provide URLs via env vars or CLI flags if defaults do not work.
   - Expected to download an abundance matrix and sample metadata.
   - Saves to data/raw/mouse_proteome/ and normalized copies to data/processed/ as needed.

2) GEO GSE225576 (paired transcriptome)
   - Downloads supplementary expression tables and sample metadata
   - Saves to data/raw/geo_gse225576/

3) PPI network for Mus musculus from STRING
   - Downloads STRING protein links for 10090 and writes a simplified edge list
     to data/ppi_edges.tsv with columns: proteinA, proteinB, score

Environment variable overrides (optional):
- MOUSE_PROTEOME_URL
- MOUSE_PROTEOME_META_URL
- GEO_GSE225576_SUPPL_URLS (comma-separated URLs)
- STRING_PPI_URL

Usage:
    python scripts/download_data.py
    python scripts/download_data.py --force
"""
from __future__ import annotations

import gzip
import io
import os
import sys
import zipfile
from pathlib import Path
from typing import Iterable, List, Optional

import requests
from tqdm import tqdm
import pandas as pd

BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
PROC_DIR = DATA_DIR / "processed"
RAW_DIR.mkdir(parents=True, exist_ok=True)
PROC_DIR.mkdir(parents=True, exist_ok=True)

# Default URLs (may change upstream; override via env if needed)
DEFAULT_MOUSE_PROTEOME_URL = os.environ.get(
    "MOUSE_PROTEOME_URL",
    # Placeholder; please set env var to the latest project URL if needed
    "https://example.org/mouse_aging_proteome/abundance_matrix.csv",
)
DEFAULT_MOUSE_PROTEOME_META_URL = os.environ.get(
    "MOUSE_PROTEOME_META_URL",
    "https://example.org/mouse_aging_proteome/sample_metadata.csv",
)
DEFAULT_STRING_PPI_URL = os.environ.get(
    "STRING_PPI_URL",
    # STRING v12 Mus musculus (10090) protein links
    "https://stringdb-static.org/download/protein.links.v12.0/10090.protein.links.v12.0.txt.gz",
)
DEFAULT_GEO_SUPPL_URLS = [
    # Override with GEO_GSE225576_SUPPL_URLS (comma-separated) if you know exact files
    # e.g., https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225576/suppl/...
]


def _download_file(url: str, dest: Path, force: bool = False) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and not force:
        return dest

    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        with tqdm(total=total, unit="B", unit_scale=True, desc=f"Downloading {dest.name}") as pbar:
            with dest.open("wb") as f:
                for chunk in r.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        f.write(chunk)
                        pbar.update(len(chunk))
    return dest


def download_mouse_proteome(force: bool = False) -> List[Path]:
    out_dir = RAW_DIR / "mouse_proteome"
    out_dir.mkdir(parents=True, exist_ok=True)

    abundance_url = DEFAULT_MOUSE_PROTEOME_URL
    metadata_url = DEFAULT_MOUSE_PROTEOME_META_URL

    abundance_path = out_dir / "abundance_matrix.csv"
    metadata_path = out_dir / "sample_metadata.csv"

    paths = []
    try:
        paths.append(_download_file(abundance_url, abundance_path, force=force))
    except Exception as e:
        print(f"Warning: failed to download abundance matrix: {e}")
    try:
        paths.append(_download_file(metadata_url, metadata_path, force=force))
    except Exception as e:
        print(f"Warning: failed to download sample metadata: {e}")

    return paths


def download_geo_gse225576(force: bool = False) -> List[Path]:
    out_dir = RAW_DIR / "geo_gse225576"
    out_dir.mkdir(parents=True, exist_ok=True)

    env_urls = os.environ.get("GEO_GSE225576_SUPPL_URLS", "").strip()
    urls: List[str] = DEFAULT_GEO_SUPPL_URLS.copy()
    if env_urls:
        urls.extend([u.strip() for u in env_urls.split(",") if u.strip()])

    if not urls:
        print("Note: No GEO supplementary URLs provided. Set GEO_GSE225576_SUPPL_URLS to download.")
        return []

    downloaded: List[Path] = []
    for url in urls:
        filename = url.split("/")[-1]
        dest = out_dir / filename
        try:
            downloaded.append(_download_file(url, dest, force=force))
            # Auto-extract simple archives
            if dest.suffix == ".zip":
                with zipfile.ZipFile(dest, "r") as zf:
                    zf.extractall(out_dir)
            elif dest.suffixes[-2:] == [".txt", ".gz"] or dest.suffix == ".gz":
                try:
                    # Leave gz as-is and let downstream handle; also write decompressed copy next to it
                    dec_path = out_dir / dest.with_suffix("").name
                    with gzip.open(dest, "rb") as f_in, dec_path.open("wb") as f_out:
                        f_out.write(f_in.read())
                except Exception:
                    pass
        except Exception as e:
            print(f"Warning: failed to download {url}: {e}")
    return downloaded


def download_string_ppi(force: bool = False) -> Path:
    out_path = DATA_DIR / "ppi_edges.tsv"
    if out_path.exists() and not force:
        return out_path

    url = DEFAULT_STRING_PPI_URL
    tmp_gz = RAW_DIR / "string_10090_links.txt.gz"
    _download_file(url, tmp_gz, force=force)

    # Parse and write simplified edge list
    # STRING format columns: protein1, protein2, combined_score
    rows = []
    with gzip.open(tmp_gz, "rt", encoding="utf-8") as fh:
        header = fh.readline().strip().split()
        col_idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            parts = line.rstrip("\n").split()
            if len(parts) < 3:
                continue
            p1 = parts[col_idx.get("protein1", 0)]
            p2 = parts[col_idx.get("protein2", 1)]
            score = parts[col_idx.get("combined_score", 2)]
            # Strip taxon prefix like "10090.ENSMUSP..." -> "ENSMUSP..."
            p1 = p1.split(".", 1)[-1] if "." in p1 else p1
            p2 = p2.split(".", 1)[-1] if "." in p2 else p2
            rows.append((p1, p2, int(score)))

    df = pd.DataFrame(rows, columns=["proteinA", "proteinB", "score"])
    df.to_csv(out_path, sep="\t", index=False)
    return out_path


def main(argv: Optional[List[str]] = None) -> None:
    import argparse

    parser = argparse.ArgumentParser(description="Download and cache datasets.")
    parser.add_argument("--force", action="store_true", help="Re-download and overwrite existing files")
    args = parser.parse_args(argv)

    print("Downloading Mouse Aging Proteome Atlas...")
    mp_paths = download_mouse_proteome(force=args.force)
    for p in mp_paths:
        print(f"  saved: {p}")

    print("\nDownloading GEO GSE225576 supplements...")
    geo_paths = download_geo_gse225576(force=args.force)
    for p in geo_paths:
        print(f"  saved: {p}")

    print("\nDownloading STRING PPI network (Mus musculus)...")
    ppi_path = download_string_ppi(force=args.force)
    print(f"  saved: {ppi_path}")

    print("\nDone. You can now run prepare_data.py")


if __name__ == "__main__":
    main()
