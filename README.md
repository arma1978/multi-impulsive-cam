# Convex Impulses Optimisation

MATLAB + C++ toolkit for impulsive collision-avoidance manoeuvre design using differential algebra (DA) propagation and convex/nonlinear optimisation.

## Overview

This project solves conjunction-avoidance problems with multiple formulations built on a shared propagation and geometry pipeline:

- SCVX (sequential convex optimisation)
- NLP (nonlinear programming)
- MILP (mixed-integer linear programming)
- Large-scale campaign analysis
- Post-processing and plotting

All calculations use km and s unless explicitly noted otherwise.

## Repository Structure

```text
convexImpulsesOptimisation/
├── main/            Top-level executable MATLAB scripts
├── src/
│   ├── routines/    Optimisation, linearisation, validation, refinement
│   ├── io/          Test-case extraction and runtime input generation
│   └── frames/      RTN/B-plane frame utilities
├── cpp/             C++ propagator sources
├── bin/             Run scripts + compiled binaries
├── utils/           MATLAB utility functions (loadPoly, evalPoly, ...)
├── input/           Curated input datasets
├── output/          Saved results and figures
├── runtime/         Regenerated scratch files (MATLAB <-> C++)
├── libdace.2.dylib
├── libdace.2.0.1.dylib
└── .gitignore
```

## Main Entry Scripts

| Script | Purpose |
|---|---|
| `main/mainSCVX.m` | Primary SCVX workflow |
| `main/mainSCVXEllipse.m` | SCVX ellipse-target variant |
| `main/mainNLP.m` | SOCP-initialised nonlinear optimisation |
| `main/mainMILP.m` | MILP formulation |
| `main/mainLargeSimSCVX.m` | Batch SCVX campaign over filtered test cases |
| `main/mainPlot.m` | Single-case plotting |
| `main/mainPlotLargeSim.m` | Large-campaign plotting |
| `main/mainBoxPlot.m` | Statistical boxplots |
| `main/mainKepJ2Compare.m` | Kepler vs J2 comparison |
| `main/mainConvert2Data.m` | Result conversion utility |
| `main/mainDataGeneration.m` | Dataset generation utility |

## Runtime Workflow

Typical execution flow:

1. Load a dataset from `input/`.
2. Select/filter a conjunction case in MATLAB (`src/io/*`).
3. Write scratch input files to `runtime/`.
4. Run a propagator binary (or `bin/run*` compile-and-run wrapper).
5. Load DA maps from `runtime/maps.dat` or `runtime/mapsRefine.dat`.
6. Build and solve the optimisation problem.
7. Save results to `output/Results*/` and figures to `output/Figures/`.

## Dependencies

### MATLAB

Validated with MATLAB R2025b on macOS.

Used toolbox/functions include:

- Optimization Toolbox (`fmincon`, `intlinprog`)
- Standard MATLAB plotting and table utilities

### MOSEK

Main optimisation scripts currently assume this MOSEK toolbox location:

```text
~/Dropbox/Work/mosek/11.1/toolbox/r2022b
```

If your environment differs, update the `addpath(...)` line near the top of `main/*.m` scripts.

### C++ Toolchain + DACE

`bin/run*` scripts compile with `g++` and link against DACE (`-ldace`).

Required dylibs are present in the repo root and in `bin/` for runtime loading on macOS.

## Data Conventions

### Persistent Inputs (`input/`)

Examples:

- `input/data.mat`
- `input/dataRev.mat`
- `input/dataEllipse.mat`
- `input/conjunctions.mat`
- `input/conjunctions.txt`
- `input/data.dat`

### Persistent Outputs (`output/`)

- `output/ResultsSCVX/`
- `output/ResultsNLP/`
- `output/ResultsMILP/`
- `output/ResultsLargeSimSCVX/`
- `output/Figures/`

### Scratch Runtime Files (`runtime/`)

Examples:

- `runtime/maps.dat`, `runtime/mapsRefine.dat`
- `runtime/xs0J2.dat`, `runtime/xs0Kep.dat`
- `runtime/xd0J2.dat`, `runtime/xd0Kep.dat`
- `runtime/xsTCA.dat`, `runtime/xdTCA.dat`
- `runtime/inputFullP.dat`, `runtime/xxs.dat`

These are regenerated and intentionally ignored by Git.

## Utility Functions (`utils/`)

- `utils/loadPoly.m`: load DA/COSY polynomial maps
- `utils/evalPoly.m`: evaluate DA polynomials at points
- `utils/po2pv.m`: orbital elements -> Cartesian state
- `utils/pv2po.m`: Cartesian state -> orbital elements
- `utils/savePlot.m`: save PNG + MATLAB `.fig`

## How To Run

Run from repository root.

### MATLAB interactive

```matlab
cd('/Users/rarm840/Dropbox/Work/CAM/convexImpulsesOptimisation')
run(fullfile('main','mainSCVX.m'))
```

### Terminal batch

```bash
cd /Users/rarm840/Dropbox/Work/CAM/convexImpulsesOptimisation
/Applications/MATLAB_R2025b.app/bin/matlab -batch "run(fullfile('main','mainSCVX.m'))"
```

Other useful commands:

```bash
/Applications/MATLAB_R2025b.app/bin/matlab -batch "run(fullfile('main','mainMILP.m'))"
/Applications/MATLAB_R2025b.app/bin/matlab -batch "run(fullfile('main','mainNLP.m'))"
/Applications/MATLAB_R2025b.app/bin/matlab -batch "quickTestMaxCases=1; run(fullfile('main','mainLargeSimSCVX.m'))"
```

For plotting with `main/mainPlot.m`, first run `main/mainSCVX.m` to generate the required result file under `output/ResultsSCVX/`.

## Rebuilding Propagators

Use wrappers in `bin/` to compile and run:

```bash
bin/runStateBackProp
bin/runPropMultiMapsFullP
bin/runPropMultiMapsKep
bin/runPropMultiMapsFullPRefine
bin/runPropMultiMapsKepRefine
```

## Validation Status

At the time of this README refresh:

- all `bin/run*` scripts were executed successfully
- all `main/*.m` scripts were executed successfully
- `mainLargeSimSCVX.m` was verified with `quickTestMaxCases=1` as a smoke test

## Notes

- Scripts assume execution from the project root.
- Some script parameters (objective selection, case index, solver settings) are intentionally top-of-file for research iteration.
- Current runtime setup is macOS-oriented because of dylib handling.

## Citation

If you use this software in academic work, please cite:

Roberto Armellin,
Collision avoidance maneuver optimization with a multiple-impulse convex formulation,
Acta Astronautica,
Volume 186,
2021,
Pages 347-362,
ISSN 0094-5765,
https://doi.org/10.1016/j.actaastro.2021.05.046.
