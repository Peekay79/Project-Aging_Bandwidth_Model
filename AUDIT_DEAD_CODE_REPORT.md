# Static Repo Usage Audit (no runtime)

Generated: 2025-10-16 06:03 UTC


Total candidates: 12  
High confidence: 2  
Medium confidence: 10

## High confidence (not reachable AND flagged by Vulture)

Top 20:
- `scripts/prepare_data.py` — LOC 247, last modified 2025-10-15 (not_reachable+vulture)
- `scripts/download_data.py` — LOC 157, last modified 2025-10-15 (not_reachable+vulture)

## Medium confidence (not reachable OR flagged by Vulture)

<details><summary>Show all</summary>


- `tools/static_audit.py` — LOC 389, last modified unknown (not_reachable)
- `scripts/run_pipeline.py` — LOC 134, last modified 2025-10-15 (not_reachable)
- `scripts/network_overlay.py` — LOC 110, last modified 2025-10-15 (not_reachable)
- `scripts/plot_headroom.py` — LOC 107, last modified 2025-10-15 (not_reachable)
- `scripts/compute_load.py` — LOC 97, last modified 2025-10-15 (not_reachable)
- `scripts/inflection.py` — LOC 96, last modified 2025-10-15 (not_reachable)
- `scripts/compute_capacity.py` — LOC 84, last modified 2025-10-15 (not_reachable)
- `scripts/compute_headroom.py` — LOC 45, last modified 2025-10-15 (not_reachable)
- `scripts/network_overlay_cli.py` — LOC 14, last modified 2025-10-15 (not_reachable)
- `scripts/inflection_cli.py` — LOC 4, last modified 2025-10-15 (not_reachable)

</details>

## Notes / Exclusions / Limitations

- Excluded directories: tests/, */venv*/, .venv/, __pycache__/, logs/, data/logs/, htmlcov/, attic/.
- Graph is built via static AST imports (import X, from X import Y) without execution.
- Imports resolving outside the repo are ignored.
- Docker compose command parsing is best-effort (regex; no YAML parser).
- Dynamic imports, reflection, and plugin loading are not detected and may cause false positives.
- Confidence: high = not reachable by import graph AND Vulture-flagged (>=80); medium = not reachable OR Vulture-flagged.
