"""PyMOL helper for coloring polymer-chain carbons by head->tail intervals.

Usage in PyMOL v2 after loading a PDB:

    run /home/joelaforet/Desktop/mupt/mupt-examples/examples_system/color_polymer_chain_carbons.py
    color_polymer_chain_carbons

Or with an explicit object name:

    color_polymer_chain_carbons ionomer_m5_production_start

Or with explicit head and tail residue names:

    color_polymer_chain_carbons output, PTH, PTT

The script defines a polymer chain as a continuous residue interval starting at
the head residue and ending at the tail residue, ordered by residue index within
each PyMOL chain ID. All atoms in those intervals are included in the generated
chain selections, but only carbon atoms are recolored. Oxygen and nitrogen retain
the standard red/blue colors.
"""

from __future__ import annotations

import colorsys
from dataclasses import dataclass

from pymol import cmd


@dataclass(frozen=True)
class PolymerInterval:
    chain: str
    start_resi: int
    end_resi: int


def _normalize_resn(resn: str) -> str:
    return resn.strip().upper()


def _safe_object_name(selection: str) -> str:
    names = cmd.get_object_list(selection)
    if not names:
        raise ValueError(f"No PyMOL object found for selection: {selection!r}")
    if len(names) > 1:
        print(f"Multiple objects matched {selection!r}; using {names[0]!r}")
    return names[0]


def _polymer_intervals(obj_name: str, head_resn: str, tail_resn: str) -> list[PolymerInterval]:
    head_resn = _normalize_resn(head_resn)
    tail_resn = _normalize_resn(tail_resn)
    residues: list[tuple[str, int, str]] = []
    seen: set[tuple[str, int]] = set()
    model = cmd.get_model(f"({obj_name})")

    for atom in model.atom:
        try:
            resi = int(atom.resi)
        except ValueError:
            continue
        key = (atom.chain, resi)
        if key in seen:
            continue
        seen.add(key)
        residues.append((atom.chain, resi, atom.resn))

    residues.sort(key=lambda item: (item[0], item[1]))

    intervals: list[PolymerInterval] = []
    active_start: tuple[str, int] | None = None
    for chain_id, resi, resn in residues:
        resn = _normalize_resn(resn)
        if resn == head_resn:
            active_start = (chain_id, resi)
        elif resn == tail_resn and active_start is not None:
            start_chain, start_resi = active_start
            if start_chain == chain_id and resi >= start_resi:
                intervals.append(PolymerInterval(chain_id, start_resi, resi))
            active_start = None

    if not intervals:
        raise ValueError(f"No {head_resn}->{tail_resn} polymer-chain intervals found")
    return intervals


def _nice_chain_rgb(index: int) -> tuple[float, float, float]:
    """Generate a pleasant non-red, non-blue chain color."""
    golden_ratio = 0.618033988749895
    hue = (0.31 + index * golden_ratio) % 1.0

    # Keep O red and N blue visually reserved.
    forbidden = (
        (0.00, 0.07),  # red
        (0.58, 0.73),  # blue/cyan-blue
        (0.93, 1.00),  # red wraparound
    )
    while any(lo <= hue <= hi for lo, hi in forbidden):
        hue = (hue + 0.11) % 1.0

    saturation = 0.58 + 0.10 * ((index * 5) % 3) / 2.0
    value = 0.78 + 0.12 * ((index * 7) % 4) / 3.0
    return colorsys.hsv_to_rgb(hue, saturation, value)


def color_polymer_chain_carbons(selection: str = "all", head_resn: str = "PEH", tail_resn: str = "PET") -> None:
    """Color carbon atoms in each head->tail polymer chain a unique color."""
    obj_name = _safe_object_name(selection)
    intervals = _polymer_intervals(obj_name, head_resn, tail_resn)

    cmd.color("red", f"({obj_name}) and elem O")
    cmd.color("blue", f"({obj_name}) and elem N")
    cmd.color("gray70", f"({obj_name}) and elem H")
    cmd.color("yelloworange", f"({obj_name}) and resn SOD")

    for idx, interval in enumerate(intervals, start=1):
        color_name = f"poly_chain_{idx:04d}"
        cmd.set_color(color_name, _nice_chain_rgb(idx - 1))
        chain_selector = f"chain {interval.chain} and " if interval.chain else ""
        selection_name = f"polymer_chain_{idx:04d}"
        chain_selection = (
            f"({obj_name}) and ({chain_selector}resi {interval.start_resi}-{interval.end_resi}) "
            "and not resn SOD"
        )
        cmd.select(selection_name, chain_selection)
        cmd.color(color_name, f"({selection_name}) and elem C")

    cmd.hide("everything", obj_name)
    cmd.show("sticks", obj_name)
    cmd.set("stick_radius", 0.12, obj_name)
    cmd.deselect()
    print(f"Colored {len(intervals)} {_normalize_resn(head_resn)}->{_normalize_resn(tail_resn)} polymer chains in {obj_name!r}")


cmd.extend("color_polymer_chain_carbons", color_polymer_chain_carbons)
