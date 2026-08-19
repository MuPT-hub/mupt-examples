"""Visualization and plotting helpers for MuPT examples."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable, Optional
import ast
import json

import numpy as np


ANGSTROM_PER_NM = 10.0


def show_available_chemistries(
    library: dict[str, dict[str, str]],
    chemistries: Optional[Iterable[str]] = None,
    max_units_per_chemistry: int = 3,
    mols_per_row: int = 3,
):
    """Draw selected repeat-unit fragments with RDKit in a notebook."""

    from IPython.display import SVG
    from rdkit import Chem
    from rdkit.Chem import Draw

    selected = list(chemistries) if chemistries is not None else sorted(library)
    mols = []
    legends = []
    skipped = []
    for chemistry in selected:
        units = library.get(chemistry, {})
        if not units:
            skipped.append(chemistry)
            continue
        for unit_name, smiles in list(units.items())[:max_units_per_chemistry]:
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if mol is None:
                skipped.append(f"{chemistry}:{unit_name}")
                continue
            mols.append(mol)
            legends.append(f"{chemistry}\n{unit_name}")
    if not mols:
        raise ValueError("No drawable repeat units were found for the selected chemistries.")
    image = Draw.MolsToGridImage(
        mols,
        legends=legends,
        molsPerRow=mols_per_row,
        subImgSize=(260, 180),
        useSVG=True,
    )
    if skipped:
        print("Skipped empty or invalid entries:", ", ".join(skipped))
    return image if isinstance(image, SVG) else SVG(image)




def _parse_metadata_value(value: Any) -> Any:
    if isinstance(value, (list, tuple, dict)):
        return value
    if not isinstance(value, str):
        return value
    for parser in (json.loads, ast.literal_eval):
        try:
            return parser(value)
        except Exception:
            pass
    return value


def box_vectors_nm_from_primitive(system: Any) -> np.ndarray:
    """Return diagonal OpenMM/OpenFF box vectors in nm from MuPT metadata."""

    params = system.metadata.get("unit_cell_parameters")
    if params is None:
        for segment in system.children:
            params = segment.metadata.get("unit_cell_parameters")
            if params is not None:
                break
    if params is None:
        raise ValueError("No unit_cell_parameters metadata found on MuPT primitive.")
    params = _parse_metadata_value(params)
    lengths_a = np.asarray(params[:3], dtype=float)
    return np.diag(lengths_a / ANGSTROM_PER_NM)


def box_lengths_a_from_primitive(system: Any) -> np.ndarray:
    """Return orthorhombic box lengths in Angstrom from MuPT metadata."""

    return np.diag(box_vectors_nm_from_primitive(system)) * ANGSTROM_PER_NM


def pdb_chain_and_resid(global_residue_idx: int) -> tuple[str, int]:
    """Return PDB-style chain/residue IDs using 9999-residue chain bins."""

    chain_ids = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    chain_idx, resid_offset = divmod(int(global_residue_idx), 9999)
    if chain_idx >= len(chain_ids):
        raise ValueError("Visualization PDBx export exceeded 26 * 9999 residues.")
    return chain_ids[chain_idx], resid_offset + 1


def rdkit_mols_from_primitive(system: Any, resname_map: dict[str, str]) -> list[Any]:
    """Export MuPT segments as RDKit molecules in role-aware atom order."""

    from mupt.interfaces.rdkit import primitive_to_rdkit_mols

    return list(primitive_to_rdkit_mols(system, resname_map=resname_map, default_atom_position=np.zeros(3)))


def make_rdkit_molecule_whole_and_centered(mol: Any, box_lengths_a: np.ndarray) -> None:
    """Unwrap one RDKit molecule and pack its centroid into the 0..L display cell."""

    from rdkit.Geometry import Point3D

    box_lengths = np.asarray(box_lengths_a, dtype=float)
    if mol.GetNumAtoms() <= 1 or box_lengths.shape != (3,) or np.any(box_lengths <= 0.0):
        return
    conf = mol.GetConformer()
    wrapped_positions = np.asarray(conf.GetPositions(), dtype=float)
    whole_positions = wrapped_positions.copy()
    adjacency = [[] for _ in range(mol.GetNumAtoms())]
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        adjacency[begin_idx].append(end_idx)
        adjacency[end_idx].append(begin_idx)

    visited = np.zeros(mol.GetNumAtoms(), dtype=bool)
    for root_idx in range(mol.GetNumAtoms()):
        if visited[root_idx]:
            continue
        visited[root_idx] = True
        stack = [root_idx]
        while stack:
            atom_idx = stack.pop()
            for neighbor_idx in adjacency[atom_idx]:
                if visited[neighbor_idx]:
                    continue
                delta = wrapped_positions[neighbor_idx] - wrapped_positions[atom_idx]
                delta -= np.round(delta / box_lengths) * box_lengths
                whole_positions[neighbor_idx] = whole_positions[atom_idx] + delta
                visited[neighbor_idx] = True
                stack.append(neighbor_idx)

    whole_positions += 0.5 * box_lengths
    centroid = whole_positions.mean(axis=0)
    whole_positions += np.mod(centroid, box_lengths) - centroid
    for atom_idx, position in enumerate(whole_positions):
        conf.SetAtomPosition(atom_idx, Point3D(float(position[0]), float(position[1]), float(position[2])))


def _atom_name(atom: Any, atom_idx: int) -> str:
    pdb_info = atom.GetPDBResidueInfo()
    if pdb_info is not None and pdb_info.GetName().strip():
        return pdb_info.GetName().strip()
    return f"{atom.GetSymbol()}{atom_idx + 1}"


def _rdkit_mols_to_openmm_topology_positions(rdkit_mols: list[Any], box_vectors_nm: np.ndarray):
    """Return an OpenMM topology and positions for visualization exports."""

    import openmm
    from openmm import unit as omm_unit
    from openmm.app import Topology, element

    topology = Topology()
    topology.setPeriodicBoxVectors(box_vectors_nm * omm_unit.nanometer)
    positions = []
    atom_lookup = {}
    residue_lookup = {}
    chain_lookup = {}
    residue_counter = 0

    for mol_idx, mol in enumerate(rdkit_mols):
        conf = mol.GetConformer()
        for atom_idx, atom in enumerate(mol.GetAtoms()):
            if atom.GetAtomicNum() == 0:
                raise ValueError("Visualization PDBx export cannot write dummy/linker atoms.")
            residue_key = (
                mol_idx,
                atom.GetProp("mupt_residue_index") if atom.HasProp("mupt_residue_index") else atom.GetProp("residue_id"),
            )
            residue = residue_lookup.get(residue_key)
            if residue is None:
                chain_id, resid = pdb_chain_and_resid(residue_counter)
                residue_counter += 1
                chain = chain_lookup.get(chain_id)
                if chain is None:
                    chain = topology.addChain(chain_id)
                    chain_lookup[chain_id] = chain
                residue_name = atom.GetProp("residue_name") if atom.HasProp("residue_name") else "UNK"
                residue = topology.addResidue(residue_name, chain, id=str(resid))
                residue_lookup[residue_key] = residue

            symbol = atom.GetSymbol()
            atom_lookup[(mol_idx, atom_idx)] = topology.addAtom(
                _atom_name(atom, atom_idx),
                element.get_by_symbol(symbol),
                residue,
            )
            position = conf.GetAtomPosition(atom_idx)
            positions.append(openmm.Vec3(position.x / ANGSTROM_PER_NM, position.y / ANGSTROM_PER_NM, position.z / ANGSTROM_PER_NM))

        for bond in mol.GetBonds():
            topology.addBond(atom_lookup[(mol_idx, bond.GetBeginAtomIdx())], atom_lookup[(mol_idx, bond.GetEndAtomIdx())])

    return topology, positions * omm_unit.nanometer


def write_rdkit_mols_to_pdb(rdkit_mols: list[Any], output_path: str | Path, box_vectors_nm: np.ndarray) -> Path:
    """Write all RDKit molecules as one visualization-friendly PDB."""

    from openmm.app import PDBFile

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    topology, positions = _rdkit_mols_to_openmm_topology_positions(rdkit_mols, box_vectors_nm)
    with output_path.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(topology, positions, handle, keepIds=True)
    return output_path


def write_rdkit_mols_to_pdbx(rdkit_mols: list[Any], output_path: str | Path, box_vectors_nm: np.ndarray) -> Path:
    """Write all RDKit molecules as one visualization-friendly PDBx/mmCIF."""

    from openmm.app import PDBxFile

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    topology, positions = _rdkit_mols_to_openmm_topology_positions(rdkit_mols, box_vectors_nm)
    with output_path.open("w", encoding="utf-8") as handle:
        PDBxFile.writeFile(topology, positions, handle)
    return output_path


def write_mupt_visualization_pdb(system: Any, resname_map: dict[str, str], output_path: str | Path) -> Path:
    """Write a centered, whole-molecule PDB structure for visualization."""

    box_lengths_a = box_lengths_a_from_primitive(system)
    rdkit_mols = rdkit_mols_from_primitive(system, resname_map=resname_map)
    for mol in rdkit_mols:
        make_rdkit_molecule_whole_and_centered(mol, box_lengths_a)
    return write_rdkit_mols_to_pdb(rdkit_mols, output_path, box_vectors_nm=box_vectors_nm_from_primitive(system))


def write_mupt_visualization_pdbx(system: Any, resname_map: dict[str, str], output_path: str | Path) -> Path:
    """Write a centered, whole-molecule mmCIF/PDBx structure for visualization."""

    box_lengths_a = box_lengths_a_from_primitive(system)
    rdkit_mols = rdkit_mols_from_primitive(system, resname_map=resname_map)
    for mol in rdkit_mols:
        make_rdkit_molecule_whole_and_centered(mol, box_lengths_a)
    return write_rdkit_mols_to_pdbx(rdkit_mols, output_path, box_vectors_nm=box_vectors_nm_from_primitive(system))




def plot_openmm_state_data(state_data_path: str | Path):
    """Load OpenMM StateDataReporter CSV and plot density and potential energy."""

    import matplotlib.pyplot as plt
    import pandas as pd

    state_df = pd.read_csv(state_data_path)
    time_column = next((col for col in state_df.columns if col.startswith('Time')), 'Step')
    density_column = next((col for col in state_df.columns if col.startswith('Density')), None)
    potential_column = next((col for col in state_df.columns if col.startswith('Potential Energy')), None)
    if density_column is None or potential_column is None:
        raise ValueError(f"State data lacks density or potential energy columns: {list(state_df.columns)}")

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.5))
    axes[0].plot(state_df[time_column], state_df[density_column], linewidth=1.5)
    axes[0].set_xlabel(time_column)
    axes[0].set_ylabel(density_column)
    axes[0].grid(alpha=0.25)
    axes[1].plot(state_df[time_column], state_df[potential_column], linewidth=1.5)
    axes[1].set_xlabel(time_column)
    axes[1].set_ylabel(potential_column)
    axes[1].grid(alpha=0.25)
    fig.tight_layout()
    return state_df, fig, axes


def plot_radius_of_gyration(rg_df: Any, y: str = "rg_nm"):
    """Plot radius of gyration versus trajectory time."""

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5.5, 3.5))
    ax.plot(rg_df["time_ps"], rg_df[y], marker="o", linewidth=1.5)
    ax.set_xlabel("time (ps)")
    ax.set_ylabel("radius of gyration (nm)" if y == "rg_nm" else "radius of gyration (A)")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    return fig, ax
