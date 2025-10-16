#!/usr/bin/env python3
from __future__ import annotations
import ast
import csv
import os
import re
import subprocess
import sys
from collections import defaultdict, deque
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

# --------- Configuration ---------
EXCLUDE_DIR_NAMES = {
    "tests",
    "__pycache__",
    "htmlcov",
    "logs",
    "attic",
}
EXCLUDE_DIR_SUBSTRINGS = {"venv", ".venv"}
EXCLUDE_DIR_PATHS = {"data/logs"}

VULTURE_EXCLUDE = "tests,venv,.venv,htmlcov,logs,data/logs,__pycache__,attic"
MIN_CONFIDENCE = 80

CSV_PATH = "dead_code_report.csv"
MD_PATH = "AUDIT_DEAD_CODE_REPORT.md"

# ---------------------------------

@dataclass(frozen=True)
class FileNode:
    path: Path  # absolute path

    @property
    def rel(self) -> Path:
        return self.path.relative_to(get_repo_root())

    @property
    def module(self) -> str:
        root = get_repo_root()
        rel = self.path.relative_to(root)
        if rel.suffix == ".py":
            rel = rel.with_suffix("")
        parts = list(rel.parts)
        return ".".join(parts)


def get_repo_root() -> Path:
    here = Path(__file__).resolve()
    # tools/ is directly under repo root
    return here.parent.parent


def should_exclude_dir(dir_path: Path) -> bool:
    rel = dir_path.relative_to(get_repo_root())
    # Fast checks by name
    name = dir_path.name
    if name in EXCLUDE_DIR_NAMES:
        return True
    for substr in EXCLUDE_DIR_SUBSTRINGS:
        if substr in name:
            return True
    # Path based exclusions
    rel_str = str(rel).replace("\\", "/")
    for ex in EXCLUDE_DIR_PATHS:
        if rel_str == ex or rel_str.startswith(ex + "/"):
            return True
    return False


def discover_python_files() -> List[Path]:
    root = get_repo_root()
    results: List[Path] = []
    for current_dir, dirnames, filenames in os.walk(root):
        # Prune excluded dirs in-place for os.walk
        pruned: List[str] = []
        for d in list(dirnames):
            full = Path(current_dir) / d
            if should_exclude_dir(full):
                pruned.append(d)
        for d in pruned:
            dirnames.remove(d)
        # Skip dot-directories like .git
        for d in list(dirnames):
            if d.startswith("."):
                dirnames.remove(d)
        for fn in filenames:
            if not fn.endswith(".py"):
                continue
            p = Path(current_dir) / fn
            # Skip files under excluded dirs (double-check)
            parts = set(p.parts)
            if any(s in p.parts for s in EXCLUDE_DIR_NAMES):
                continue
            # Skip dot-directories
            if any(part.startswith(".") for part in p.parts):
                continue
            results.append(p.resolve())
    return results


def build_module_index(py_files: Iterable[Path]) -> Dict[str, Set[Path]]:
    """Map module dotted names -> file paths that provide them.

    Strategy:
    - For every file a/b/c.py => module "a.b.c"
    - For every package a/b/__init__.py => module "a.b"
    Both map to sets to avoid collisions; prefer leaf modules during resolution.
    """
    root = get_repo_root()
    index: Dict[str, Set[Path]] = defaultdict(set)
    for f in py_files:
        rel = f.relative_to(root)
        if rel.name == "__init__.py":
            pkg = ".".join(rel.parent.parts)
            if pkg:
                index[pkg].add(f)
        else:
            mod = ".".join(rel.with_suffix("").parts)
            if mod:
                index[mod].add(f)
    return index


def parse_imports(file_path: Path) -> List[Tuple[str, Optional[int]]]:
    """Parse a file and extract import targets as dotted names.

    Returns list of (import_target, preference_rank) where lower rank is preferred.
    preference_rank is used to mildly prefer explicit module files over package __init__ fallbacks.
    """
    try:
        text = file_path.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return []
    try:
        tree = ast.parse(text, filename=str(file_path))
    except Exception:
        return []

    root = get_repo_root()
    current_rel = file_path.relative_to(root)
    current_mod_parts = list(current_rel.with_suffix("").parts)

    imports: List[Tuple[str, Optional[int]]] = []

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                name = alias.name  # e.g., "pkg.mod"
                if name:
                    imports.append((name, None))
        elif isinstance(node, ast.ImportFrom):
            level = getattr(node, "level", 0) or 0
            module = node.module  # may be None for "from . import x"

            # Build base module parts considering relative level
            if level > 0:
                base_parts = current_mod_parts[:-level]
            else:
                base_parts = []
            if module:
                base_parts = base_parts + module.split(".")

            for alias in node.names:
                target = alias.name  # may be a module or a symbol
                # Prefer module target if it exists; otherwise, try just the base module
                dotted_candidate = ".".join(base_parts + ([target] if target else []))
                if dotted_candidate:
                    imports.append((dotted_candidate, 0))
                # Also consider the base module (package) itself as a light fallback
                if base_parts:
                    imports.append((".".join(base_parts), 1))
    return imports


def resolve_import_to_files(dotted: str, module_index: Dict[str, Set[Path]]) -> List[Path]:
    """Resolve a dotted module to repository files.

    Resolution tries exact match. We do not attempt parent shortening beyond package fallback
    (handled by parse_imports). Returns files only if present in the repo.
    """
    if dotted in module_index:
        # Prefer leaf module files (non-__init__) by sorting
        files = list(module_index[dotted])
        files.sort(key=lambda p: 0 if p.name != "__init__.py" else 1)
        return files
    return []


def build_import_graph(py_files: List[Path]) -> Tuple[Dict[Path, Set[Path]], Dict[str, Set[Path]]]:
    module_index = build_module_index(py_files)
    edges: Dict[Path, Set[Path]] = defaultdict(set)

    for f in py_files:
        imports = parse_imports(f)
        # Prioritize lower preference rank first, and dedupe
        seen: Set[str] = set()
        for dotted, _rank in sorted(imports, key=lambda x: (x[1] is not None, x[1] or 0)):
            if dotted in seen:
                continue
            seen.add(dotted)
            targets = resolve_import_to_files(dotted, module_index)
            for t in targets:
                if t != f:
                    edges[f].add(t)
    # Ensure all nodes exist in edges dict
    for f in py_files:
        edges.setdefault(f, set())
    return edges, module_index


def discover_entrypoints(py_files: List[Path]) -> Set[Path]:
    root = get_repo_root()
    candidates: Set[Path] = set()

    # 1) startup_discord_interface.py
    sdi = root / "startup_discord_interface.py"
    if sdi.exists():
        candidates.add(sdi.resolve())

    # 2) pods/memory/app.py or any pods/*/app.py
    pods_dir = root / "pods"
    if pods_dir.exists() and pods_dir.is_dir():
        mem_app = pods_dir / "memory" / "app.py"
        if mem_app.exists():
            candidates.add(mem_app.resolve())
        for child in pods_dir.rglob("app.py"):
            candidates.add(child.resolve())

    # 3) docker-compose*.yml command entries pointing to python files (best effort regex)
    for yml in root.glob("docker-compose*.yml"):
        try:
            txt = yml.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            continue
        # Look for lines like: command: python -m pkg.mod or command: python script.py
        for m in re.finditer(r"command:\s*(.+)", txt):
            cmd = m.group(1)
            # simple split by whitespace
            parts = re.split(r"\s+", cmd.strip())
            # Direct script invocation
            for p in parts:
                if p.endswith(".py"):
                    script = (yml.parent / p).resolve() if not os.path.isabs(p) else Path(p)
                    if script.exists():
                        candidates.add(script.resolve())
            # python -m module
            if "-m" in parts:
                try:
                    idx = parts.index("-m")
                    mod = parts[idx + 1]
                    # Try mapping module to files via module index (approx): construct expected path
                    mod_path = root.joinpath(*mod.split(".")).with_suffix(".py")
                    if mod_path.exists():
                        candidates.add(mod_path.resolve())
                    init_path = root.joinpath(*mod.split(".")) / "__init__.py"
                    if init_path.exists():
                        candidates.add(init_path.resolve())
                except Exception:
                    pass

    # Filter to only files present in py_files
    allowed = {p.resolve() for p in py_files}
    return {p for p in candidates if p in allowed}


def reachable_from(entrypoints: Set[Path], edges: Dict[Path, Set[Path]]) -> Set[Path]:
    if not entrypoints:
        return set()
    visited: Set[Path] = set()
    dq: deque[Path] = deque()
    for ep in entrypoints:
        if ep in edges:
            dq.append(ep)
            visited.add(ep)
    while dq:
        cur = dq.popleft()
        for nxt in edges.get(cur, ()): 
            if nxt not in visited:
                visited.add(nxt)
                dq.append(nxt)
    return visited


def run_vulture(repo_root: Path) -> Set[Path]:
    cmd = [
        sys.executable,
        "-m",
        "vulture",
        ".",
        "--min-confidence",
        str(MIN_CONFIDENCE),
        "--exclude",
        VULTURE_EXCLUDE,
    ]
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(repo_root),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="ignore",
        )
    except FileNotFoundError:
        return set()

    out = proc.stdout or ""
    flagged: Set[Path] = set()
    # Typical lines: path:line: message (NN% confidence)
    for line in out.splitlines():
        line = line.strip()
        if not line:
            continue
        # vulture also prints summary lines; ignore those that have no ':'
        if ":" not in line:
            continue
        path_part = line.split(":", 1)[0]
        # Ignore non-py
        if not path_part.endswith(".py"):
            continue
        p = (repo_root / path_part).resolve()
        if p.exists():
            flagged.add(p)
    return flagged


def count_loc(file_path: Path) -> int:
    try:
        with file_path.open("r", encoding="utf-8", errors="ignore") as f:
            count = 0
            for raw in f:
                s = raw.strip()
                if not s:
                    continue
                if s.startswith("#"):
                    continue
                count += 1
            return count
    except Exception:
        return 0


def git_last_modified(file_path: Path) -> str:
    try:
        proc = subprocess.run(
            ["git", "log", "-1", "--format=%cs", "--", str(file_path)],
            cwd=str(get_repo_root()),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="ignore",
            timeout=10,
        )
        out = (proc.stdout or "").strip()
        return out or "unknown"
    except Exception:
        return "unknown"


def write_csv(rows: List[Dict[str, str | int]], path: Path) -> None:
    fieldnames = ["path", "confidence", "reason", "loc", "last_modified"]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


def write_markdown(
    repo_root: Path,
    high: List[Dict[str, str | int]],
    medium: List[Dict[str, str | int]],
    path: Path,
) -> None:
    now = datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")
    total = len(high) + len(medium)
    top20 = high[:20]

    def fmt_row(r: Dict[str, str | int]) -> str:
        return f"- `{r['path']}` — LOC {r['loc']}, last modified {r['last_modified']} ({r['reason']})"

    lines: List[str] = []
    lines.append(f"# Static Repo Usage Audit (no runtime)\n")
    lines.append(f"Generated: {now}\n")
    lines.append("")
    lines.append(f"Total candidates: {total}  ")
    lines.append(f"High confidence: {len(high)}  ")
    lines.append(f"Medium confidence: {len(medium)}\n")

    lines.append("## High confidence (not reachable AND flagged by Vulture)\n")
    if high:
        lines.append("Top 20:")
        lines.extend(fmt_row(r) for r in top20)
    else:
        lines.append("(none)")

    lines.append("")
    lines.append("## Medium confidence (not reachable OR flagged by Vulture)\n")
    if medium:
        lines.append("<details><summary>Show all</summary>\n")
        lines.append("")
        for r in medium:
            lines.append(fmt_row(r))
        lines.append("")
        lines.append("</details>\n")
    else:
        lines.append("(none)\n")

    lines.append("## Notes / Exclusions / Limitations\n")
    lines.append("- Excluded directories: tests/, */venv*/, .venv/, __pycache__/, logs/, data/logs/, htmlcov/, attic/.")
    lines.append("- Graph is built via static AST imports (import X, from X import Y) without execution.")
    lines.append("- Imports resolving outside the repo are ignored.")
    lines.append("- Docker compose command parsing is best-effort (regex; no YAML parser).")
    lines.append("- Dynamic imports, reflection, and plugin loading are not detected and may cause false positives.")
    lines.append("- Confidence: high = not reachable by import graph AND Vulture-flagged (>=80); medium = not reachable OR Vulture-flagged.")

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    repo_root = get_repo_root()
    py_files = discover_python_files()

    edges, module_index = build_import_graph(py_files)
    entrypoints = discover_entrypoints(py_files)
    reachable = reachable_from(entrypoints, edges)

    vulture_flagged = run_vulture(repo_root)

    # Prepare candidate rows
    def make_row(f: Path) -> Optional[Dict[str, str | int]]:
        is_reachable = f in reachable
        is_flagged = f in vulture_flagged
        if not is_reachable or is_flagged:
            reason_parts = []
            if not is_reachable:
                reason_parts.append("not_reachable")
            if is_flagged:
                reason_parts.append("vulture")
            reason = "+".join(reason_parts)
            confidence = "high" if (not is_reachable and is_flagged) else "medium"
            rel = f.relative_to(repo_root)
            return {
                "path": str(rel).replace("\\", "/"),
                "confidence": confidence,
                "reason": reason,
                "loc": count_loc(f),
                "last_modified": git_last_modified(f),
            }
        return None

    rows: List[Dict[str, str | int]] = []
    for f in py_files:
        r = make_row(f)
        if r is not None:
            rows.append(r)

    # Sort: high first, then LOC desc, then path
    rows.sort(key=lambda r: (0 if r["confidence"] == "high" else 1, -(r["loc"] or 0), str(r["path"])) )

    # Split for markdown sections
    high_rows = [r for r in rows if r["confidence"] == "high"]
    medium_rows = [r for r in rows if r["confidence"] == "medium"]

    write_csv(rows, repo_root / CSV_PATH)
    write_markdown(repo_root, high_rows, medium_rows, repo_root / MD_PATH)

    print(f"Wrote {CSV_PATH} with {len(rows)} rows")
    print(f"Wrote {MD_PATH}")
    if entrypoints:
        print("Entrypoints detected:")
        for ep in sorted(entrypoints):
            print(f" - {ep.relative_to(repo_root)}")
    else:
        print("No entrypoints detected.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
