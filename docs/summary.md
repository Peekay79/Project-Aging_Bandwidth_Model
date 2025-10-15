# Headroom analysis summary

## Data sources and versions

- Mouse Aging Proteome Atlas (TMT proteomics) — downloaded via scripts/download_data.py (see docs/ASSUMPTIONS.md for exact files).
- GEO GSE225576 transcriptome (optional; cached in data/raw/geo_gse225576/).
- STRING PPI (Mus musculus, v12) — simplified to data/ppi_edges.tsv.

## Normalisation steps

- Within-batch median normalization across TMT channels (per batch_id).
- Per-tissue z-score for each gene across samples (ages).
- Missingness filter: retain genes quantified in ≥70% of samples per tissue.

## Capacity/Load formulas

- Capacity per sample = z-scored sum of normalized abundances for chaperone ∪ UPS ∪ autophagy genes (from data/protein_annotations_seed.csv).
- Load per sample = z-scored sum of aggregation-prone and stress-proxy genes; fallback includes top 10% by disorder/missingness proxy when explicit markers are absent.
- Headroom = Capacity − Load (z-space).

## Headline observations

- Brain: headroom decreases with age (slope=-0.045).
- Heart: headroom decreases with age (slope=-0.038).
- Intestine: headroom decreases with age (slope=-0.013).
- Kidney: headroom decreases with age (slope=-0.042).
- Liver: headroom decreases with age (slope=-0.107).
- Lung: headroom decreases with age (slope=-0.080).
- Muscle: headroom increases with age (slope=0.102).
- Spleen: headroom decreases with age (slope=-0.101).

Earliest headroom < 0 observed in Muscle at ~3 months.

Earliest estimated inflections (AIC-selected):
- Spleen: logistic inflection at ~2.6 months
- Liver: logistic inflection at ~2.8 months
- Muscle: logistic inflection at ~3.2 months
- Brain: piecewise inflection at ~18.7 months
- Kidney: piecewise inflection at ~19.1 months

## Limitations and next steps

- Real proteome URLs may change; configure via env vars. ID mapping is best-effort using mygene.
- Load fallback uses disorder/missingness proxy; refine with curated stress markers when available.
- Consider batch-specific adjustments beyond median normalization and integrate transcriptome pairing.
