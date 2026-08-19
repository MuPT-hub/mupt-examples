"""Workshop helpers for dense polymer melt MuPT demos.

These functions intentionally live in the examples repository while the public
MuPT API is still evolving. They keep the notebooks focused on the workflow:
choose chemistry, build a dense melt with AA-DPD, hand off to OpenMM, and analyze
the trajectory with MDAnalysis.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Optional
import ast
import json

import numpy as np

from recipes.saamr import MeltBuildPlan


ANGSTROM_PER_NM = 10.0


@dataclass(frozen=True)
class MDRunConfig:
    """OpenMM settings in MD-user terms instead of OpenMM object terms."""

    temperature_k: float = 450.0
    pressure_atm: float = 1.0
    timestep_fs: float = 1.0
    friction_per_ps: float = 1.0
    nvt_time_ns: float = 0.01
    npt_time_ns: float = 0.02
    frames_to_save: int = 25
    force_field: str = "openff-2.2.1.offxml"
    charge_method: str = "openff-gnn-am1bcc-1.0.0.pt"
    platform_name: Optional[str] = None
    platform_properties: dict[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class OpenMMRunResult:
    """Paths and energies from the minimal OpenMM handoff workflow."""

    trajectory_path: Path
    nvt_trajectory_path: Path
    npt_trajectory_path: Path
    nvt_state_data_path: Path
    npt_state_data_path: Path
    final_pdb_path: Path
    initial_energy: Any
    minimized_energy: Any
    final_energy: Any
    nvt_steps: int
    npt_steps: int


def find_examples_root(start: Optional[Path] = None) -> Path:
    """Locate the mupt-examples repository root from a notebook cwd."""

    start = Path.cwd() if start is None else Path(start)
    for candidate in (start, *start.parents):
        if (candidate / "examples_system").exists() and (candidate / "README.md").exists():
            return candidate
    raise RuntimeError("Could not locate the mupt-examples repository root.")


def load_repeat_unit_libraries(paths: Optional[Iterable[Path]] = None) -> dict[str, dict[str, str]]:
    """Load bundled repeat-unit SMILES libraries."""

    root = find_examples_root()
    if paths is None:
        paths = [
            root / "examples_system" / "repeat_unit_SMILES.json",
            root / "examples_system" / "repeat_units_arom_hetcycles.json",
        ]
    library: dict[str, dict[str, str]] = {}
    for path in paths:
        with Path(path).open("r", encoding="utf-8") as handle:
            data = json.load(handle)
        for chemistry, units in data.items():
            if chemistry in library and library[chemistry] and units:
                merged = dict(library[chemistry])
                merged.update(units)
                library[chemistry] = merged
            else:
                library[chemistry] = units
    return library


def show_available_chemistries(
    library: dict[str, dict[str, str]],
    chemistries: Optional[Iterable[str]] = None,
    max_units_per_chemistry: int = 3,
    mols_per_row: int = 3,
):
    """Draw selected repeat-unit fragments with RDKit in a notebook."""

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
    image = Draw.MolsToGridImage(mols, legends=legends, molsPerRow=mols_per_row, subImgSize=(260, 180))
    if skipped:
        print("Skipped empty or invalid entries:", ", ".join(skipped))
    return image


def summarize_build_plan(plan: MeltBuildPlan):
    """Return a compact pandas summary table for a planned melt."""

    import pandas as pd

    return pd.DataFrame(
        [
            ("chemistry", plan.recipe_name),
            ("sizing mode", plan.sizing_mode),
            ("chains", plan.n_chains),
            ("repeat units", int(sum(plan.chain_lengths))),
            ("chain length range", f"{min(plan.chain_lengths)}-{max(plan.chain_lengths)}"),
            ("box lengths (nm)", tuple(round(x, 3) for x in plan.box_lengths_nm)),
            ("target density (g/cm^3)", round(plan.target_density_g_cm3, 4)),
            ("actual density (g/cm^3)", round(plan.actual_density_g_cm3, 4)),
            ("actual mass (Da)", round(plan.actual_mass_da, 2)),
        ],
        columns=["quantity", "value"],
    )


def summarize_aa_dpd(system: Any, result: Any):
    """Return a compact pandas summary of an AA-DPD run."""

    import pandas as pd

    summary = system.metadata.get("all_atom_dpd_summary", {})
    return pd.DataFrame(
        [
            ("converged", result.converged),
            ("steps", result.steps),
            ("elapsed (s)", round(result.elapsed_s, 2)),
            ("atoms", summary.get("n_atoms", len(result.atoms))),
            ("bonds", summary.get("n_bonds", len(result.bonds))),
            ("angles", summary.get("n_angles", len(result.angles))),
            ("dihedrals", summary.get("n_dihedrals", len(result.dihedrals))),
            ("box lengths (A)", tuple(round(x, 3) for x in result.box_lengths_a)),
        ],
        columns=["quantity", "value"],
    )


def write_demo_mupt_sdf(system: Any, path: str | Path, resname_map: dict[str, str]) -> Path:
    """Write a role-aware temporary `.mupt.sdf` and return its actual path."""

    from mupt.temporary.sdf import MUPT_SDF_SUFFIX, write_primitive_to_sdf

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    write_primitive_to_sdf(system, path, resname_map=resname_map, default_atom_position=np.zeros(3))
    if str(path).endswith(MUPT_SDF_SUFFIX):
        return path
    if path.suffix == ".sdf":
        return path.with_suffix(MUPT_SDF_SUFFIX)
    return Path(f"{path}{MUPT_SDF_SUFFIX}")


def load_demo_mupt_sdf(path: str | Path, sanitize: bool = False) -> Any:
    """Load a temporary `.mupt.sdf` into a role-aware MuPT hierarchy."""

    from mupt.temporary.sdf import primitive_from_mupt_sdf

    return primitive_from_mupt_sdf(path, sanitize=sanitize)


def write_manifest(path: str | Path, **entries: Any) -> Path:
    """Write a small JSON manifest connecting the two workshop notebooks."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    def clean(value: Any) -> Any:
        if isinstance(value, Path):
            return str(value)
        if isinstance(value, (np.integer, np.floating)):
            return value.item()
        if isinstance(value, tuple):
            return list(value)
        if isinstance(value, list):
            return [clean(v) for v in value]
        if isinstance(value, dict):
            return {str(k): clean(v) for k, v in value.items()}
        return value

    path.write_text(json.dumps(clean(entries), indent=2) + "\n", encoding="utf-8")
    return path


def load_manifest(path: str | Path) -> dict[str, Any]:
    """Load a workshop manifest JSON file."""

    return json.loads(Path(path).read_text(encoding="utf-8"))


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


def openff_molecules_from_primitive(system: Any, resname_map: dict[str, str]) -> list[Any]:
    """Convert MuPT segments to OpenFF Molecule objects with coordinates."""

    from openff.toolkit import Molecule

    return [
        Molecule.from_rdkit(mol, allow_undefined_stereo=True, hydrogens_are_explicit=True)
        for mol in rdkit_mols_from_primitive(system, resname_map=resname_map)
    ]


def unique_openff_molecules(molecules: list[Any]) -> list[Any]:
    """Return isomorphically unique OpenFF Molecules for preset charges."""

    return list(dict.fromkeys(molecules))


def assign_openff_charges(molecules: list[Any], charge_method: str) -> None:
    """Assign charges in-place, using NAGL explicitly for `.pt` models."""

    if charge_method in {"zeros", "zero"}:
        from openff.units import unit

        for molecule in molecules:
            molecule.partial_charges = np.zeros(molecule.n_atoms) * unit.elementary_charge
        return
    if charge_method.endswith(".pt") or charge_method.startswith("openff-gnn"):
        from openff.toolkit.utils import ToolkitRegistry
        from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper

        if not NAGLToolkitWrapper.is_available():
            raise RuntimeError(
                "The default charge method requires openff-nagl. Install the "
                "mupt-demo-env environment or choose charge_method='zeros' only "
                "for quick debugging."
            )
        registry = ToolkitRegistry([NAGLToolkitWrapper()])
        for molecule in molecules:
            molecule.assign_partial_charges(charge_method, toolkit_registry=registry)
        return
    for molecule in molecules:
        molecule.assign_partial_charges(charge_method)


def build_openff_interchange(
    system: Any,
    resname_map: dict[str, str],
    force_field: str = "openff-2.2.1.offxml",
    charge_method: str = "openff-gnn-am1bcc-1.0.0.pt",
) -> Any:
    """Build an OpenFF Interchange from a MuPT primitive."""

    from openff.toolkit import ForceField, Topology
    from openff.units import unit

    molecules = openff_molecules_from_primitive(system, resname_map=resname_map)
    charge_molecules = unique_openff_molecules(molecules)
    assign_openff_charges(charge_molecules, charge_method)
    topology = Topology.from_molecules(molecules)
    positions_a = np.vstack([mol.conformers[0].m_as(unit.angstrom) for mol in molecules])
    topology.box_vectors = box_vectors_nm_from_primitive(system) * unit.nanometer
    interchange = ForceField(force_field).create_interchange(topology, charge_from_molecules=charge_molecules)
    interchange.positions = positions_a * unit.angstrom
    interchange.box = topology.box_vectors
    return interchange


def _interchange_to_openmm_system(interchange: Any) -> Any:
    for method_name in ("to_openmm_system", "to_openmm"):
        method = getattr(interchange, method_name, None)
        if method is None:
            continue
        try:
            return method(combine_nonbonded_forces=True)
        except TypeError:
            return method()
    raise RuntimeError("OpenFF Interchange cannot be converted to an OpenMM System.")


def _interchange_to_openmm_topology(interchange: Any) -> Any:
    method = getattr(interchange, "to_openmm_topology", None)
    if method is not None:
        return method()
    return interchange.topology.to_openmm()


def _steps(time_ns: float, timestep_fs: float) -> int:
    return int(round(time_ns * 1_000_000.0 / timestep_fs))


def run_openmm_workflow(interchange: Any, config: MDRunConfig, output_dir: str | Path, prefix: str) -> OpenMMRunResult:
    """Run minimization, NVT, and NPT with OpenMM and save separate trajectories."""

    import openmm
    from openmm import LangevinMiddleIntegrator, MonteCarloBarostat, unit as omm_unit
    from openmm.app import DCDReporter, PDBFile, Simulation, StateDataReporter

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    nvt_trajectory_path = output_dir / f"{prefix}_openmm_nvt.dcd"
    npt_trajectory_path = output_dir / f"{prefix}_openmm_npt.dcd"
    nvt_state_data_path = output_dir / f"{prefix}_nvt_state_data.csv"
    npt_state_data_path = output_dir / f"{prefix}_npt_state_data.csv"
    final_pdb_path = output_dir / f"{prefix}_final.pdb"

    system = _interchange_to_openmm_system(interchange)
    topology = _interchange_to_openmm_topology(interchange)
    if getattr(interchange, "box", None) is not None and hasattr(topology, "setPeriodicBoxVectors"):
        topology.setPeriodicBoxVectors(interchange.box.to_openmm())

    integrator = LangevinMiddleIntegrator(
        config.temperature_k * omm_unit.kelvin,
        config.friction_per_ps / omm_unit.picosecond,
        config.timestep_fs * omm_unit.femtosecond,
    )
    if config.platform_name:
        platform = openmm.Platform.getPlatformByName(config.platform_name)
        simulation = Simulation(topology, system, integrator, platform, config.platform_properties)
    else:
        simulation = Simulation(topology, system, integrator)
    print(f"OpenMM platform: {simulation.context.getPlatform().getName()}")
    simulation.context.setPositions(interchange.positions.to_openmm())
    if getattr(interchange, "box", None) is not None:
        simulation.context.setPeriodicBoxVectors(*interchange.box.to_openmm())

    initial_state = simulation.context.getState(getEnergy=True)
    initial_energy = initial_state.getPotentialEnergy()
    simulation.minimizeEnergy()
    minimized_state = simulation.context.getState(getEnergy=True)
    minimized_energy = minimized_state.getPotentialEnergy()

    nvt_steps = _steps(config.nvt_time_ns, config.timestep_fs)
    npt_steps = _steps(config.npt_time_ns, config.timestep_fs)
    def close_reporter(reporter: Any) -> None:
        close = getattr(getattr(reporter, "_out", None), "close", None)
        if close is not None:
            close()

    state_reporter_kwargs = dict(
        step=True,
        time=True,
        potentialEnergy=True,
        kineticEnergy=True,
        temperature=True,
        volume=True,
        density=True,
        speed=True,
    )

    if nvt_steps:
        nvt_interval = max(1, nvt_steps // max(1, config.frames_to_save))
        nvt_reporter = DCDReporter(str(nvt_trajectory_path), nvt_interval)
        nvt_state_reporter = StateDataReporter(
            str(nvt_state_data_path),
            nvt_interval,
            **state_reporter_kwargs,
        )
        simulation.reporters.extend([nvt_reporter, nvt_state_reporter])
        simulation.context.setVelocitiesToTemperature(config.temperature_k * omm_unit.kelvin)
        simulation.step(nvt_steps)
        for reporter in (nvt_reporter, nvt_state_reporter):
            simulation.reporters.remove(reporter)
            close_reporter(reporter)
    if npt_steps:
        system.addForce(
            MonteCarloBarostat(
                config.pressure_atm * omm_unit.atmosphere,
                config.temperature_k * omm_unit.kelvin,
                25,
            )
        )
        simulation.context.reinitialize(preserveState=True)
        npt_interval = max(1, npt_steps // max(1, config.frames_to_save))
        npt_reporter = DCDReporter(str(npt_trajectory_path), npt_interval)
        npt_state_reporter = StateDataReporter(
            str(npt_state_data_path),
            npt_interval,
            **state_reporter_kwargs,
        )
        simulation.reporters.extend([npt_reporter, npt_state_reporter])
        simulation.step(npt_steps)
        for reporter in (npt_reporter, npt_state_reporter):
            simulation.reporters.remove(reporter)
            close_reporter(reporter)

    final_state = simulation.context.getState(getEnergy=True, getPositions=True, enforcePeriodicBox=False)
    final_energy = final_state.getPotentialEnergy()
    with final_pdb_path.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(topology, final_state.getPositions(), handle)
    return OpenMMRunResult(
        trajectory_path=npt_trajectory_path,
        nvt_trajectory_path=nvt_trajectory_path,
        npt_trajectory_path=npt_trajectory_path,
        nvt_state_data_path=nvt_state_data_path,
        npt_state_data_path=npt_state_data_path,
        final_pdb_path=final_pdb_path,
        initial_energy=initial_energy,
        minimized_energy=minimized_energy,
        final_energy=final_energy,
        nvt_steps=nvt_steps,
        npt_steps=npt_steps,
    )


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
