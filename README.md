# Plasma Proteomics Replicability — R Sample (.rda)

This repository provides **minimal, runnable R scripts** to demonstrate association replicability analysis (Directional Consistency, DC) using your `.rda` files.

## Data
Your uploaded files have been placed in `data/`:
- `test.rda`
- `train.rda`
- `Z_Example.rda`

> Note: Do **not** commit platform-specific binaries (e.g., `*.dll`). They are ignored via `.gitignore`.

## Quick start (R)
```r
# 0) Inspect .rda contents (object names / structure)
source("R/00_inspect_rda.R")

# 1) Compute DC (edit file names/columns inside to match your objects)
source("R/01_dc_replicability.R")
```

## What you need to edit
Open `R/01_dc_replicability.R` and set:
- Which `.rda` files to load (defaults to `train.rda` and `test.rda`).
- The object names inside each `.rda` (e.g., `assocA`, `assocB`).
- Column names used for merge and signs (default assumes: `protein`, `phenotype`, `beta`).

## Acknowledgement
If you drew inspiration from the structure of the proteome–phenome atlas, please add a link to their GitHub in your published README.

## License
This code is released under the MIT License. See the `LICENSE` file for details.