# Legacy notebooks (MuPT 0.1.3a1)

These notebooks use MuPT's pre-refactor API and the old built-in AA-DPD builder. Run them in [`conda-envs/legacy.yml`](../conda-envs/legacy.yml).
They are kept unchanged here while the team decides, for each one, whether to port it, fold it into a workflow, or retire it.

| Notebook | What it shows | Notes |
|---|---|---|
| `1_Building_Dense_Polymer_Melts` | dense melts from monomer SMILES, DOP, sequence statistics and density (old AA-DPD builder + `recipes/saamr.py`) | builder being retired (MuPT-hub/mupt#139) |
| `2_Simulating_Dense_Polymer_Melts` | OpenFF parameterization, export, short OpenMM run, MDAnalysis | reads notebook 1's outputs |
| `3_Analysis_Examples_Dense_Polymer_Melts` | analysis of the melts | reads notebook 1's outputs |
| `5_CG_Systems` | coarse-grained workflows through flowerMD | source for `workflows/03_paam_cg_dpd` |
| `mupt_to_flower-PAAm` | CG bead-spring DPD for polyacrylamide | source for `workflows/03_paam_cg_dpd` |
| `mupt_to_flower-dendr` | CG PAMAM dendrimer | needs a dendrimer generator on the refactored API |
| `hierarchy_on_peptides` | peptide hierarchy from a FASTA string | currently pulled into the MuPT docs as a tutorial |
| `multiscale_chain_placement` | concurrent all-atom / CG chain placement | |
| `mdfiles_with_openff` | force field parameterization and MD file export | |
