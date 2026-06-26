# mupt-examples

## Installation
To obtain a copy of these examples, clone this repo and navigate into it, via:
```sh
git clone https://github.com/MuPT-hub/mupt-examples # this repo
cd mupt-examples
```
Then set up a virtual environment with the Multiscale Polymer Toolkit installed to run the notebooks bundled here. The [`mupt` toolkit repo](https://github.com/MuPT-hub/mupt) provides details on how to do so, but if you don't like clicking you can also run the commands below, assuming you have a package management systems such as [Mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) (recommended) or [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) installed on your machine:
```sh
mamba env create -f conda-envs/release-env.yml
mamba activate mupt-env
git clone https://github.com/MuPT-hub/mupt # the toolkit 
cd mupt
pip install .
```
You can then safely remove the `mupt` directory from your file system if you wish to

## Index of examples

### End-to-End workflow examples (quickstart):
0. [Building a dense polymer melt](./examples_system/Building_Dense_Polymer_Melts.ipynb)
Allows a user to provide various monomer SMILES, DOP, statistical distribution of monomers, and a density to construct initial coordinates for a dense polymer melt.
1. [Simulate the melt](./examples_system/Simulating_Dense_Polymer_Melts.ipynb)
Showcases how one can take the system constructed in 0. and assign force field parameters via OpenFF Interchange, export to OpenMM|GROMACS|LAMMPS, and demo a short OpenMM simulation, with optional analysis via MDAnalysis.

### System building examples (high-level):
0. [Simple peptide example](./examples_system/hierarchy_on_peptides.ipynb)  
Assembles a MuPT representation component hierarchy for a peptide generated from your choice of FASTA string  
Depicts the topological and hierarchical organization of parts in a colorful graphic with 2D coordinates
1. [Backmapped chain/melt building](./examples_system/ellipsoidal_chain_placement.ipynb)  
Shows how to perform a backmapping-like build procedure to go from nothing to an MD starting structure using MuPT tools  
Allows specification of arbitrary homopolymer or copolymer chemistry, number of chains, and distribution of chain length  
2. [Parameterization, export, and MD simulations with OpenFF](./examples_system/mdfiles_with_openff.ipynb)  
Shows how to use the OpenFF toolkit to parameterize and export MD input files from a MuPT-build starting structure (obtained from the previous demo)  
Concludes by running a short MD simulation in OpenMM

### MuPT Representation components tutorials (low-level):
0. ["Molecule from scratch" - bare basics](./examples_repr/mol_from_scratch_basic.ipynb)  
Shows how to initialize a basic "bond graph" hierarchy with no coordinates
1. ["Molecule from scratch" - advanced usage](./examples_repr/mol_from_scratch_advanced.ipynb)  
"Kitchen sink" tutorial of how to build a molecule fragment with complete hierarchical info using nothing but MuPT representation components  
Generates physically-plausible coordinates and topology using no force fields or external coordinate embedding algorithms (just basic knowledge of point groups)  
and shows how one can leverage the concurrent, arbitrary resolutions structure of the representation to organize initialization around chemically-meaningful motifs   
2. [Molecule from SMILES/RDKit](./examples_repr/mupt_primitives_rdkit.ipynb)  
Shows how to instantiate MuPT representation hierarchies from RDKit Mols and SMILES string
3. [Rigid body shapes](./examples_repr/mupt_shapes.ipynb)  
Shows off that different kinds of shapes that MuPT Primitive objects can take on (PointClouds, Spheres, Ellipsoids, and Rods)  
4. [Inter-body spatial alignment](./examples_repr/connector_alignment.ipynb)
Generates animated trajectories which illustrate the different kinds of inter-body alignment operations MuPT uses to spatially align all MuPT representation components

<!-- ## DPD notebooks:
The DPD notebooks are meant to demonstrate a system initialization workflow of polymer various systems and configuration routines. These configurations are relaxed using a HOOMD simulation with DPD forces. A list of the notebooks and their descriptions:
+ DPD_init.ipynb: initializes a polymer system using a random walk and relaxes the system using DPD forces.
+ DPD_init_ell.ipynb: initializes a random configuration of ellipsoid chains, based on the flowerMD ellipsoid model, and relaxes the system using DPD forces.
+ DPD_init_mbuild.ipynb: initializes a polymer system using an mbuild2.0 random walk and relaxes the system using DPD forces.

Environment for  DPD_init.ipynb can be installed with `dpd-env.yml` or with:

```
conda create -n mupt-init -c conda-forge matplotlib scipy gsd freud hoomd fresnel jupyter
```

As of Aug 1 2025 DPD_init_ell.ipynb requires a developer install of [flowerMD](https://github.com/cmelab/flowermd)

As of Sep 25 2025 DPD_init_mbuild.ipynb requires a developer install of [mbuild](https://github.com/mosdef-hub/mbuild) `gh pr checkout 1261` -->

