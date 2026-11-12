# mupt-examples

Example notebooks for [MuPT](https://github.com/MuPT-hub/mupt), the multiscale polymer toolkit.

> These notebooks use MuPT's **refactored API**, which is still in review on the [MuPT PR stack](https://github.com/MuPT-hub/mupt/pulls) (#141 through #166). Environments pin the top branch of that stack until it is released. Notebooks for the previous release (0.1.3a1) are in [`legacy/`](./legacy).

## Installation

Clone the examples, then set up an environment with **pixi** (fastest) or **conda/mamba**. Both give the same packages.

```sh
git clone https://github.com/MuPT-hub/mupt-examples
cd mupt-examples
```

### Option A: pixi

[pixi](https://pixi.sh) creates the environment inside the repo and pins exact versions in `pixi.lock`, so everyone runs the same thing. The pixi environments currently support **Linux** (and Windows through [WSL](https://learn.microsoft.com/windows/wsl/install)); on macOS, use Option B for now. Install pixi once:

```sh
curl -fsSL https://pixi.sh/install.sh | sh    # then open a new terminal
```

Then start JupyterLab. The first run creates the environment (a few minutes); later runs start right away:

```sh
pixi run lab
```

That's it: notebooks open with every package available. To use one of the smaller environments instead, add `-e <name>`:

```sh
pixi run -e workflows lab
```

**Using VS Code or a JupyterLab you already have?** Register the environment as a Jupyter kernel once, then pick "MuPT examples (default)" from the kernel menu:

```sh
pixi run install-kernel                # or: pixi run -e workflows install-kernel
```

### Option B: conda / mamba

```sh
mamba env create -f conda-envs/mupt-examples.yml
mamba activate mupt-examples
jupyter lab
```

To use the environment from VS Code or another Jupyter install: `python -m ipykernel install --user --name mupt-examples` (with the environment active).

### Environments

| pixi (`-e <name>`) | conda file | Use it for | Includes |
|---|---|---|---|
| `default` | [`mupt-examples.yml`](./conda-envs/mupt-examples.yml) | everything except `legacy/` | all of the below |
| `representation` | [`representation.yml`](./conda-envs/representation.yml) | `representation/`, `interfaces/` | MuPT, RDKit, plotting |
| `workflows` | [`workflows.yml`](./conda-envs/workflows.yml) | `workflows/` | + flowerMD, HOOMD, mBuild/GMSO, OpenFF, OpenMM, MDAnalysis |
| `workflows-hoomd5`, `workflows-hoomd5-gpu-cuda11`, `workflows-hoomd5-gpu-cuda12` | | `workflows/` on HOOMD 5.2 (see below) | as `workflows` |
| `legacy` | [`legacy.yml`](./conda-envs/legacy.yml) | `legacy/` | MuPT 0.1.3a1 (pre-refactor API) |

### Older GPUs: HOOMD 5.2

The default environments use the latest HOOMD. If your GPU or NVIDIA driver is too old for it, use HOOMD 5.2. Check which CUDA version your driver supports with `nvidia-smi` (top right, "CUDA Version"), then:

| Your driver supports | Run |
|---|---|
| CUDA 11.8 – 12.5 | `pixi run -e workflows-hoomd5-gpu-cuda11 lab` |
| CUDA 12.6 or newer | `pixi run -e workflows-hoomd5-gpu-cuda12 lab` |
| no NVIDIA GPU | `pixi run -e workflows-hoomd5 lab` (CPU build) |

GPU builds are Linux-only; WSL supports NVIDIA GPUs, so Windows users can run them there.

## Index

Status is against the refactored API. "Docs" is the proposed [Diátaxis](https://diataxis.fr) category for the MuPT documentation site ([MuPT-hub/mupt#105](https://github.com/MuPT-hub/mupt/issues/105)).

### Workflows: end to end, from SMILES to MD input files

| Notebook | What it shows | Status | Docs |
|---|---|---|---|
| [`01_pe_all_atom_dpd`](./workflows/01_pe_all_atom_dpd.ipynb) | polyethylene melt: MuPT → mBuild → flowerMD AA-DPD → back onto MuPT → PDB + GROMACS/LAMMPS/OpenMM inputs (OpenFF Interchange) | ✅ runs | Tutorial |
| [`02_ps_all_atom_openmm`](./workflows/02_ps_all_atom_openmm.ipynb) | isotactic polystyrene melt (10 chains × 20): same pipeline with stereocenters kept through AA-DPD, then OpenFF Sage + a short OpenMM run | ✅ runs | Tutorial |
| [`03_paam_cg_dpd`](./workflows/03_paam_cg_dpd.ipynb) | coarse-grained polyacrylamide (5 chains × 20 beads): one bead per repeat unit, flowerMD bead-spring DPD, back onto MuPT, CG PDB + HOOMD gsd | ✅ runs | Tutorial |

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
pixi run test-workflows                       # or, with conda: cd workflows && bash notebook_tests.sh
pixi run -e workflows-hoomd5 test-workflows   # same notebooks on HOOMD 5.2
```
