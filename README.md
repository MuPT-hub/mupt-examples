# mupt-examples

Example notebooks for [MuPT](https://github.com/MuPT-hub/mupt), the multiscale polymer toolkit.

> These notebooks use MuPT's **refactored API**, which is still in review on the [MuPT PR stack](https://github.com/MuPT-hub/mupt/pulls) (#141 through #166). Environments pin the top branch of that stack until it is released. Notebooks for the previous release (0.1.3a1) are in [`legacy/`](./legacy).

## Installation

```sh
git clone https://github.com/MuPT-hub/mupt-examples
cd mupt-examples
mamba env create -f conda-envs/mupt-examples.yml
mamba activate mupt-examples
jupyter lab
```

| Environment | Use it for | Includes |
|---|---|---|
| [`mupt-examples.yml`](./conda-envs/mupt-examples.yml) | everything except `legacy/` | all of the below |
| [`representation.yml`](./conda-envs/representation.yml) | `representation/`, `interfaces/` | MuPT, RDKit, plotting |
| [`workflows.yml`](./conda-envs/workflows.yml) | `workflows/` | + flowerMD, HOOMD, mBuild/GMSO, OpenFF, OpenMM, MDAnalysis |
| [`legacy.yml`](./conda-envs/legacy.yml) | `legacy/` | MuPT 0.1.3a1 (pre-refactor API) |

## Index

Status is against the refactored API. "Docs" is the proposed [Diátaxis](https://diataxis.fr) category for the MuPT documentation site ([MuPT-hub/mupt#105](https://github.com/MuPT-hub/mupt/issues/105)).

### Workflows: end to end, from SMILES to MD input files

| Notebook | What it shows | Status | Docs |
|---|---|---|---|
| [`01_pe_all_atom_dpd`](./workflows/01_pe_all_atom_dpd.ipynb) | polyethylene melt: MuPT → mBuild → flowerMD AA-DPD → back onto MuPT → PDB + GROMACS/LAMMPS/OpenMM inputs (OpenFF Interchange) | ✅ runs | Tutorial |
| [`02_ps_all_atom_openmm`](./workflows/02_ps_all_atom_openmm.ipynb) | isotactic polystyrene melt (10 chains × 20): same pipeline with stereocenters kept through AA-DPD, then OpenFF Sage + a short OpenMM run | ✅ runs | Tutorial |
| `03_paam_cg_dpd` | coarse-grained polyacrylamide bead-spring DPD via flowerMD | planned | Tutorial |

### Representation: how MuPT describes molecules

| Notebook | What it shows | Status | Docs |
|---|---|---|---|
| [`01_hierarchy_basics`](./representation/01_hierarchy_basics.ipynb) | building a bond-graph hierarchy with no coordinates | needs port (`visualize_topology` removed) | Tutorial |
| [`02_hierarchy_advanced`](./representation/02_hierarchy_advanced.ipynb) | "kitchen sink": a fragment with full hierarchy and coordinates from point groups | needs port; likely replaced by simpler examples | Not in docs |
| [`03_shapes`](./representation/03_shapes.ipynb) | rigid-body shapes: point clouds, spheres, ellipsoids, rods | imports OK, not re-run | How-To |
| [`04_connector_alignment`](./representation/04_connector_alignment.ipynb) | animations of the alignment operations between bodies | imports OK, not re-run | How-To |

### Interfaces: getting molecules in and out of MuPT

| Notebook | What it shows | Status | Docs |
|---|---|---|---|
| [`01_rdkit`](./interfaces/01_rdkit.ipynb) | MuPT hierarchies from RDKit molecules and SMILES | needs port (`visualize_topology` removed) | How-To |

### Analysis

Planned: analysis examples moving over from [mupt-analysis](https://github.com/MuPT-hub/mupt-analysis) ([#22](https://github.com/MuPT-hub/mupt-examples/issues/22)).

### Legacy (MuPT 0.1.3a1)

Kept unchanged while we decide where each one belongs; see [`legacy/README.md`](./legacy/README.md). Not part of the docs site for now.

## Testing

`notebook_tests_master.sh` runs the notebooks in `representation/`, `interfaces/` and `legacy/` (MuPT's CI runs this script). Workflow notebooks need flowerMD and HOOMD, so run them separately:

```sh
cd workflows && bash notebook_tests.sh
```
