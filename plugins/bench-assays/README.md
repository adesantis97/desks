# bench-assays

Quantify bench assays from photos, with a human in the loop.

## Skills

- `colony-count` — colony formation assay (crystal violet)
  - Finds the dishes and the colonies, numbers them on an overlay
  - Review page: click to remove or add a colony, save, recompute
  - Output: `counts.csv` (colonies, mean size, fold change vs CTRL), `colonies.csv` (one row per colony)
- `invasion-count` — Transwell invasion / migration (crystal violet)
  - Measures the % of each field covered by stain, not a cell count
  - Review page: one threshold slider per condition, click to exclude bad fields
  - Output: per field, per condition, and a summary with fold change and p-value over replicates

## Requirements

- Python 3.9+ with `requirements.txt`
- R + tidyverse only for the optional figures (`--plot`)

## Validation

- Synthetic data (`examples/make_examples.py`): colony counts exact on 2-dish and 6-dish layouts; invasion cell counts exact on 24 synthetic fields (fold change 0.412 vs 0.412). Real membranes: confluent fields (> 60% area) give a minimum count.
- Real data: built and tuned on melanoma and glioblastoma cell lines, crystal-violet stained, phone and microscope photos.
- Not validated against a blinded manual count: always look at the overlay before trusting a number.

## Folder layout

```
skills/colony-count/     SKILL.md · colony_count.py · review_template.html · plot_colonies.R
skills/invasion-count/   SKILL.md · invasion_count.py · review_template.html · plot_invasion.R
examples/                make_examples.py (synthetic test images with ground truth)
```
