"""OpenFF and OpenMM workflow helpers for MuPT examples."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import numpy as np

from utilities.visualization import box_vectors_nm_from_primitive


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




def openff_molecules_from_primitive(system: Any, resname_map: dict[str, str]) -> list[Any]:
    """Convert MuPT segments to OpenFF Molecule objects with coordinates."""

    from openff.toolkit import Molecule

    from mupt.interfaces.rdkit import primitive_to_rdkit_mols

    rdkit_mols = primitive_to_rdkit_mols(
        system,
        resname_map=resname_map,
        default_atom_position=np.zeros(3),
    )
    return [
        Molecule.from_rdkit(mol, allow_undefined_stereo=True, hydrogens_are_explicit=True)
        for mol in rdkit_mols
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



