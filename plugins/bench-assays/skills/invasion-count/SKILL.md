---
name: invasion-count
description: |
  Measure an invasion/migration assay (Transwell, stained cells photographed under the microscope)
  from several fields per condition, with an interactive review: Claude proposes the threshold, the
  user moves it and drops the bad fields, Claude recomputes.
  Primary measure = % of the field covered by the stain (invaded area), not the number of cells:
  on a confluent membrane single cells cannot be separated.
  Outputs fields.csv (one value per photo), condition.csv (mean ± SD of the fields) and, over the
  whole tree, invasion_summary.csv with fold change vs CTRL across biological replicates.
  TRIGGER: photos of stained Transwell membranes; "invasion assay"; "migration assay"; `/invasion-count`.
  NOT for: colonies on dishes (use `/colony-count`); crystal violet OD read in wells.
---

# invasion-count

`SKILL_DIR` = the folder that contains this SKILL.md. Dependencies: `pip install numpy pandas pillow scipy scikit-image` (R + tidyverse only for figures).

## Where it runs
- **Claude Code** (the user's own machine): everything as written below. `review` opens the page in their browser.
- **Cowork** (or any sandbox that cannot open the user's browser): work inside the folder the user shared, so they can see the files.
  - Install the dependencies first; if `pip install` refuses, add `--user` or `--break-system-packages`.
  - Do **not** run `review` (the page would be served inside the sandbox, out of the user's reach). Instead tell the user to open `<ROOT>/<rep>/<line>/<cond>/_invasion/review.html` from their folder with a double click.
  - On that page "Save" copies the corrections to the clipboard: the user pastes them in the chat, you write them to `<ROOT>/<rep>/<line>/<cond>/_invasion/corrections.json` and continue as after a save.
  - No R in the sandbox → skip figures unless `Rscript` exists.

## The folder tree
```
<ROOT>/<replicate>/<line>/<condition>/*.bmp   ← the 20X photos, ~10 per condition
                                     _invasion/  ← written by the skill, inside the condition
```
Names are free: only the three levels matter. Photos at another magnification (`..._10X.bmp`) are
skipped automatically: use `--mag 10` if those are the ones you need. Files without a `NNX` tag are always used.

## Steps
1. Put the photo tree in the project folder.
2. Ask the user only what the names do not tell: which condition is the CTRL, which lines to compare.
3. Measure everything:
   ```
   python3 "$SKILL_DIR/invasion_count.py" run <ROOT>
   ```
4. **Open the review of two conditions first** (one dense, one sparse: calibrate on the pair):
   ```
   python3 "$SKILL_DIR/invasion_count.py" review <ROOT>/<rep>/<line>/<cond>
   ```
   The page: threshold slider (green = what counts as cell), click a field to exclude it,
   mean ± SD and CV at the top. "Save and recompute" writes `_invasion/corrections.json`.
5. Take the threshold the user approved and **re-run everything with it**, no more pages:
   ```
   python3 "$SKILL_DIR/invasion_count.py" run <ROOT> --thr <threshold>
   ```
6. Aggregate and show the table:
   ```
   python3 "$SKILL_DIR/invasion_count.py" aggregate <ROOT> --ctrl CTRL
   ```
7. Figures only if the user asks: `--plot [--colour <colour>]`.
8. Look at 2-3 `_invasion/*_overlay.jpg` yourself before saying the number is good.

## The threshold is the only choice
- Within the same replicate **and** line the threshold must be **the same for every condition**:
  that is what makes the percentages comparable. Default 0.10, calibrated on 20X crystal-violet membranes.
- Across replicates it may change (different staining and exposure): the fold change vs the CTRL of
  the same replicate absorbs the difference. That is why aggregation is done on fold changes, not percentages.
- A wrong threshold shows in the overlay: green on the beige membrane = too low, holes inside cells = too high.

## Files
- `_invasion/fields.csv` — one row per photo: `area_pct`, `iod`, `used`, `threshold`.
- `_invasion/condition.csv` — mean, SD and CV% of the kept fields.
- `_invasion/corrections.json` — `{"threshold": 0.12, "exclude": ["photo.bmp", ...]}`.
- `_invasion/<photo>_overlay.jpg` — the mask in green over the photo. Always look at it.
- `<ROOT>/invasion_fields.csv · invasion_conditions.csv · invasion_summary.csv` — from `aggregate`.
- `area_pct` = % of the field above threshold, after removing objects < 20 px (debris, pores).
- `iod` = mean stain index over the field: independent of the threshold, a cross-check.
- `fold_change` = mean area of the condition / mean area of the CTRL **of the same replicate and line**.
- `p_vs_ctrl` = one-sample t-test of the fold changes against 1, across biological replicates (needs n ≥ 3).

## Known limits
- Confluent field (> 60% area): the measure saturates and differences shrink. Flag it.
- Out-of-focus or very dark cells turn brownish and the G−B index catches them poorly: check the overlay.
- Air bubbles, folds and tears in the membrane: exclude them by hand on the page, do not correct them.
- Different number of fields per condition (9 vs 10): the mean holds, but write it in the report.
- It is not a cell count: never write "n cells", write "invaded area (%)".
