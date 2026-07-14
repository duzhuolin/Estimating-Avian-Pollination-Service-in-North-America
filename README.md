# Code and data for: "Uneven avian population trends compromise resilience in pollination service across North America"

This repository contains all data, models, and scripts required to reproduce the analyses, figures, and tables reported in the manuscript.

## Repository structure

```
├── Main script.R                                          # Primary analysis script (all figures & tables)
├── input/
│   ├── avian nectarivore data for computing FVR.csv       # Species traits for pollination capacity
│   ├── Rosenberg et al species list.csv                   # Species metadata from Rosenberg et al. (2019)
│   └── Rosenberg et al annual indices of abundance.csv    # BBS annual abundance indices
├── models/
│   ├── Rosenberg et al model.txt                          # Hierarchical model (main)
│   ├── GAM model smoothing indices jagam.txt              # GAM smoothing of abundance indices
│   ├── tempgam.txt / tempgampred.txt                      # GAM diagnostics
│   └── uncertainty source analysis-*.txt                  # One-at-a-time uncertainty partitioning models (7 files)
├── renv/                                                  # R environment (renv)
│   └── renv.lock                                          # Package lockfile for exact reproducibility
├── Estimating Avian Pollination Service in North America.Rproj
├── LICENSE                                                # GPL-3.0
└── README.md
```

## Data

### Species traits (`input/avian nectarivore data for computing FVR.csv`)

Trait data for 17 North American nectarivorous bird species, including body mass, allometric coefficients for field metabolic rate (FMR), diet fractions (FN, FNE), nectar traits (NV, NC), service days, degree of ecological specialization, and migratory behavior. See Supplementary Data 1 in the manuscript for variable definitions and data sources.

### Population data (`input/Rosenberg et al *.csv`)

Annual abundance indices and population estimates from the North American Breeding Bird Survey (BBS), processed following Rosenberg et al. (2019, *Science*). The species list is subset to the 17 nectarivores analyzed in this study.

## Models

Bayesian hierarchical models implemented in JAGS, adapted from Rosenberg et al. (2019). The main model (`Rosenberg et al model.txt`) estimates annual population size for each species from BBS indices. The uncertainty-partitioning models run the full Monte Carlo pipeline with only one variable allowed to vary at a time, quantifying each input variable's contribution to total uncertainty in continental AFV change.