# plasma-proteomics-replicability

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22952620.svg)](https://doi.org/10.5281/zenodo.22952620)

This repository provides **minimal, runnable R scripts** for demonstrating the association replicability analysis based on Directional Consistency (DC) described in the study:

**"Merits and challenges of plasma proteomics on association replicability"**

A fixed release of the demonstration code, **v1.0.0**, is permanently archived on Zenodo:  https://doi.org/10.5281/zenodo.22952620

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

The code in this repository is released under the MIT License. See the `LICENSE` file for details.

## Citation

If you use this code, please cite the archived version:

Jiao, Z., Zhang, Y., Lai, Y., Kang, J., Ma, L., Zhao, W., You, J., Cheng, W. & Feng, J. *Merits and challenges of plasma  proteomics on association replicability*. plasma-proteomics-replicability, Zenodo, https://doi.org/10.5281/zenodo.22952620 (2026).
