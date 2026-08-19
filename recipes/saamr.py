"""Reusable SAAMR recipes for statistical linear polymer melts."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any
import json

import networkx as nx
import numpy as np


DA_PER_NM3_TO_G_CM3 = 1.0 / 602.214076
ANGSTROM_PER_NM = 10.0


@dataclass(frozen=True)
class StatisticalLinearPolymerRecipe:
    """User-selected repeat-unit chemistry and chain composition."""

    name: str
    repeat_unit_smiles: dict[str, str]
    head_unit: str
    tail_unit: str
    mid_unit_distribution: dict[str, float]
    resname_map: dict[str, str]


@dataclass(frozen=True)
class MeltBuildPlan:
    """Concrete system-size plan for a dense melt build."""

    recipe_name: str
    sizing_mode: str
    target_density_g_cm3: float
    box_lengths_nm: tuple[float, float, float]
    chain_lengths: list[int]
    chain_sequences: list[list[str]]
    target_mass_da: float
    actual_mass_da: float
    actual_density_g_cm3: float
    random_seed: int

    @property
    def box_lengths_a(self) -> tuple[float, float, float]:
        return tuple(ANGSTROM_PER_NM * length for length in self.box_lengths_nm)

    @property
    def n_chains(self) -> int:
        return len(self.chain_lengths)


def validate_recipe(recipe: StatisticalLinearPolymerRecipe) -> None:
    """Check that a recipe references known repeat units and probabilities."""

    known = set(recipe.repeat_unit_smiles)
    requested = {recipe.head_unit, recipe.tail_unit, *recipe.mid_unit_distribution}
    missing = sorted(requested - known)
    if missing:
        raise ValueError(f"Recipe references repeat units not present in repeat_unit_smiles: {missing}")
    total = sum(recipe.mid_unit_distribution.values())
    if not np.isclose(total, 1.0):
        raise ValueError(f"Middle repeat-unit probabilities must sum to 1.0; got {total}.")
    for unit_name in requested:
        if unit_name not in recipe.resname_map:
            raise ValueError(f"Recipe resname_map lacks a three-letter code for {unit_name!r}.")


def build_repeat_unit_lexicon(recipe: StatisticalLinearPolymerRecipe, semiminor_fraction: float = 0.5) -> dict[str, Any]:
    """Build oriented MuPT residue templates from recipe SMILES."""

    validate_recipe(recipe)
    from mupt.geometry.coordinates.reference import CoordAxis, origin
    from mupt.geometry.shapes import Ellipsoid
    from mupt.geometry.transforms.rigid import rigid_vector_coalignment
    from mupt.interfaces.rdkit import suppress_rdkit_logs
    from mupt.interfaces.smiles import primitive_from_smiles
    from mupt.roles import PrimitiveRole

    lexicon = {}
    axis = CoordAxis.X
    with suppress_rdkit_logs():
        for unit_name, smiles in recipe.repeat_unit_smiles.items():
            if unit_name not in {recipe.head_unit, recipe.tail_unit, *recipe.mid_unit_distribution}:
                continue
            unit = primitive_from_smiles(
                smiles,
                ensure_explicit_Hs=True,
                embed_positions=True,
                label=unit_name,
            )
            unit.role = PrimitiveRole.RESIDUE
            unit.metadata["residue_name"] = recipe.resname_map[unit_name]
            for atom in unit.children:
                atom.role = PrimitiveRole.PARTICLE

            mapped_atoms = unit.search_hierarchy_by(
                lambda prim: "molAtomMapNumber" in prim.metadata,
                min_count=2,
            )
            head_atom, tail_atom = mapped_atoms[:2]
            head_pos = head_atom.shape.centroid
            tail_pos = tail_atom.shape.centroid
            major_radius = np.linalg.norm(tail_pos - head_pos) / 2.0
            axis_vec = np.zeros(3, dtype=float)
            axis_vec[axis.value] = major_radius
            unit.rigidly_transform(
                rigid_vector_coalignment(
                    vector1_start=head_pos,
                    vector1_end=tail_pos,
                    vector2_start=origin(3),
                    vector2_end=axis_vec,
                    t1=0.5,
                    t2=0.0,
                )
            )
            radii = np.full(3, semiminor_fraction * major_radius)
            radii[axis.value] = major_radius
            unit.shape = Ellipsoid(radii)
            lexicon[unit_name] = unit
    return lexicon


def primitive_mass_da(primitive: Any) -> float:
    """Return total atom mass for a MuPT primitive."""

    mass = 0.0
    for atom in primitive.leaves:
        if atom.element is None:
            raise ValueError(f"Atom primitive {atom.label!r} has no element.")
        mass += float(atom.element.mass)
    return mass


def sequence_repeat_units(recipe: StatisticalLinearPolymerRecipe, chain_length: int, rng: np.random.Generator) -> list[str]:
    """Sample one chain sequence from a recipe."""

    if chain_length < 2:
        raise ValueError("chain_length must be at least 2 for head and tail units.")
    mid_names = list(recipe.mid_unit_distribution)
    probabilities = list(recipe.mid_unit_distribution.values())
    mids = rng.choice(mid_names, size=chain_length - 2, p=probabilities).astype(object).tolist()
    return [recipe.head_unit, *mids, recipe.tail_unit]


def chain_mass_da(sequence: list[str], lexicon: dict[str, Any]) -> float:
    """Return exact atom mass for one chain sequence."""

    return sum(primitive_mass_da(lexicon[unit_name]) for unit_name in sequence)


def expected_repeat_unit_mass_da(recipe: StatisticalLinearPolymerRecipe, lexicon: dict[str, Any]) -> float:
    """Return composition-weighted middle repeat-unit mass for box-fill planning."""

    return sum(
        float(probability) * primitive_mass_da(lexicon[unit_name])
        for unit_name, probability in recipe.mid_unit_distribution.items()
    )


def _sample_chain_lengths(
    n_chains: int,
    chain_length_range: tuple[int, int],
    rng: np.random.Generator,
) -> list[int]:
    low, high = chain_length_range
    if n_chains < 1:
        raise ValueError("n_chains must be positive.")
    if low < 2 or high < low:
        raise ValueError("chain_length_range must satisfy 2 <= min <= max.")
    return rng.integers(low, high + 1, size=n_chains).astype(int).tolist()


def _density_g_cm3(mass_da: float, box_lengths_nm: tuple[float, float, float]) -> float:
    volume_nm3 = float(np.prod(box_lengths_nm))
    return mass_da * DA_PER_NM3_TO_G_CM3 / volume_nm3


def _box_for_density_nm(mass_da: float, density_g_cm3: float) -> tuple[float, float, float]:
    if density_g_cm3 <= 0.0:
        raise ValueError("density_g_cm3 must be positive.")
    volume_nm3 = mass_da * DA_PER_NM3_TO_G_CM3 / density_g_cm3
    length_nm = float(volume_nm3 ** (1.0 / 3.0))
    return (length_nm, length_nm, length_nm)


def plan_box_from_density(
    recipe: StatisticalLinearPolymerRecipe,
    lexicon: dict[str, Any],
    n_chains: int,
    chain_length_range: tuple[int, int],
    density_g_cm3: float,
    random_seed: int = 42,
) -> MeltBuildPlan:
    """Workflow A: choose chemistry, chain count, and density; compute box."""

    rng = np.random.default_rng(random_seed)
    chain_lengths = _sample_chain_lengths(n_chains, chain_length_range, rng)
    sequences = [sequence_repeat_units(recipe, length, rng) for length in chain_lengths]
    actual_mass = sum(chain_mass_da(sequence, lexicon) for sequence in sequences)
    box_lengths_nm = _box_for_density_nm(actual_mass, density_g_cm3)
    return MeltBuildPlan(
        recipe_name=recipe.name,
        sizing_mode="density_to_box",
        target_density_g_cm3=float(density_g_cm3),
        box_lengths_nm=box_lengths_nm,
        chain_lengths=chain_lengths,
        chain_sequences=sequences,
        target_mass_da=actual_mass,
        actual_mass_da=actual_mass,
        actual_density_g_cm3=_density_g_cm3(actual_mass, box_lengths_nm),
        random_seed=random_seed,
    )


def plan_chains_for_box(
    recipe: StatisticalLinearPolymerRecipe,
    lexicon: dict[str, Any],
    box_lengths_nm: tuple[float, float, float],
    density_g_cm3: float,
    chain_length_range: tuple[int, int],
    random_seed: int = 42,
) -> MeltBuildPlan:
    """Workflow B: choose chemistry, density, and box; compute chain lengths."""

    from mupt.builders.all_atom_dpd import AllAtomDPDBuilder

    repeat_mass = expected_repeat_unit_mass_da(recipe, lexicon)
    box_lengths_a = tuple(ANGSTROM_PER_NM * length for length in box_lengths_nm)
    plan = AllAtomDPDBuilder.plan_uniform_chain_lengths_for_box(
        density_g_cm3=density_g_cm3,
        box_lengths_a=box_lengths_a,
        repeat_unit_mass_amu=repeat_mass,
        chain_length_min=chain_length_range[0],
        chain_length_max=chain_length_range[1],
        random_seed=random_seed,
    )
    rng = np.random.default_rng(random_seed)
    sequences = [sequence_repeat_units(recipe, length, rng) for length in plan.chain_lengths]
    actual_mass = sum(chain_mass_da(sequence, lexicon) for sequence in sequences)
    return MeltBuildPlan(
        recipe_name=recipe.name,
        sizing_mode="box_density_to_chain_count",
        target_density_g_cm3=float(density_g_cm3),
        box_lengths_nm=tuple(float(length) for length in box_lengths_nm),
        chain_lengths=list(plan.chain_lengths),
        chain_sequences=sequences,
        target_mass_da=float(plan.target_mass_amu),
        actual_mass_da=actual_mass,
        actual_density_g_cm3=_density_g_cm3(actual_mass, box_lengths_nm),
        random_seed=random_seed,
    )


def build_statistical_linear_polymer_melt(recipe: StatisticalLinearPolymerRecipe, lexicon: dict[str, Any], plan: MeltBuildPlan) -> Any:
    """Build a role-aware MuPT hierarchy from a concrete melt plan."""

    from mupt.mupr.primitives import Primitive
    from mupt.mupr.topology import TopologicalStructure
    from mupt.roles import PrimitiveRole

    universe = Primitive(label=f"{recipe.name}_dense_melt", role=PrimitiveRole.UNIVERSE)
    universe.metadata.update(
        {
            "recipe_name": recipe.name,
            "target_density_g_cm3": str(plan.target_density_g_cm3),
            "actual_density_g_cm3": str(plan.actual_density_g_cm3),
            "box_lengths_nm": json.dumps(list(plan.box_lengths_nm)),
            "unit_cell_parameters": json.dumps([*plan.box_lengths_a, 90.0, 90.0, 90.0]),
        }
    )
    for chain_idx, sequence in enumerate(plan.chain_sequences):
        segment = Primitive(label=f"chain_{chain_idx:04d}", role=PrimitiveRole.SEGMENT)
        for repeat_idx, unit_name in enumerate(sequence):
            residue = lexicon[unit_name].copy()
            residue.label = f"{unit_name}_{repeat_idx:03d}"
            residue.role = PrimitiveRole.RESIDUE
            residue.metadata.update(
                {
                    "residue_name": recipe.resname_map[unit_name],
                    "chain_index": str(chain_idx),
                    "repeat_index": str(repeat_idx),
                }
            )
            for atom in residue.children:
                atom.role = PrimitiveRole.PARTICLE
            segment.attach_child(residue)
        segment.set_topology(
            nx.path_graph(segment.children_by_handle.keys(), create_using=TopologicalStructure),
            max_registration_iter=100,
        )
        universe.attach_child(segment)
    return universe
