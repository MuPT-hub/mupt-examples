"""Shared notebook plumbing for MuPT examples."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable, Optional
import json

import numpy as np

from recipes.saamr import MeltBuildPlan


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



