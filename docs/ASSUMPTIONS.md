# Assumptions and processing notes

This document is appended by `scripts/prepare_data.py` when run. Key points:

- Input files are auto-parsed with delimiter and header detection. If both abundance and metadata are missing, a mock-to-tidy fallback is used purely for smoke tests.
- Protein IDs are mapped to gene symbols using a provided gene column when present; otherwise via `mygene` best-effort mapping and a local cache in `data/processed/idmap_cache.csv`.
- Normalisation: within-batch median normalisation across TMT channels (per `batch_id`), then per-tissue z-scoring for each gene across samples.
- Missingness filter: retain genes quantified in at least 70% of samples per tissue.
- Capacity per sample = z-scored sum of chaperone ∪ UPS ∪ autophagy genes.
- Load per sample = z-scored sum of aggregation-prone and stress-proxy genes; fallback includes top 10% of genes by disorder/missingness proxy.
- Headroom = Capacity − Load (z-space). Aggregation to tissue×age uses the mean with 95% CI.
- Batch IDs discovered during preparation are recorded below when available.
## Data preparation assumptions

- Protein identifier column inferred as `protein_id`.
- No explicit gene symbol column found; used ID mapping to gene symbols.
- Within-batch median normalization applied across TMT channels (per `batch_id`).
- Detected batch IDs: batch_1.
- Per-tissue z-scoring applied for each gene across samples (ages).
- Missingness filter: kept genes quantified in ≥70% of samples per tissue.
- Metadata columns used: `sample_id, tissue, age_months, batch_id`.
- ID mapping unresolved for 31978 proteins (example: P0001, P0002, P0003, P0004, P0005, P0006, P0007, P0008, P0009, P0010…). Fallback: keep original ID as `gene_symbol`.

## Data preparation assumptions

- Protein identifier column inferred as `protein_id`.
- No explicit gene symbol column found; used ID mapping to gene symbols.
- Within-batch median normalization applied across TMT channels (per `batch_id`).
- Detected batch IDs: batch_1.
- Per-tissue z-scoring applied for each gene across samples (ages).
- Missingness filter: kept genes quantified in ≥70% of samples per tissue.
- Metadata columns used: `sample_id, tissue, age_months, batch_id`.

