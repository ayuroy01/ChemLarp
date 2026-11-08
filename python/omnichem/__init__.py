"""Python bindings for OmniChem Core via the C ABI (include/omnichem/omnichem_c.h).

No compiled extension is needed: the shared library is loaded with ctypes.
Set OMNICHEM_LIB to override the library path.

    >>> import omnichem
    >>> omnichem.parse("CCO")["structure"]["species"][0]["formula"]
    'C2H6O'
    >>> sys = omnichem.System()
    >>> sys.add("CCO")
    0
    >>> sys.atom_x()        # zero-copy numpy view into the SoA column
"""

from __future__ import annotations

import ctypes
import json
import os
import platform
from pathlib import Path
from typing import Sequence

__all__ = ["parse", "canonical", "kinetics", "check_reaction", "equilibrate", "System", "OmniChemError", "library_path"]


class OmniChemError(RuntimeError):
    def __init__(self, diagnostics: list[dict]):
        self.diagnostics = diagnostics
        lines = [f"{d['code']} @ {d.get('position')}: {d['message']}" for d in diagnostics if d["severity"] == "error"]
        super().__init__("\n".join(lines) or "OmniChem error")


def _candidate_libraries() -> list[Path]:
    build = Path(__file__).resolve().parents[2] / "build"
    system = platform.system()
    if system == "Windows":
        names = ["omnichem.dll", "libomnichem.dll"]
    elif system == "Darwin":
        names = ["libomnichem.dylib"]
    else:
        names = ["libomnichem.so"]
    # CMake multi-config generators (Visual Studio, Xcode) add a per-config subfolder.
    dirs = [build, build / "Release", build / "RelWithDebInfo", build / "Debug"]
    return [d / n for d in dirs for n in names]


def library_path() -> Path:
    """Shared library to load: $OMNICHEM_LIB, else the first build output that exists."""
    env = os.environ.get("OMNICHEM_LIB")
    if env:
        return Path(env)
    candidates = _candidate_libraries()
    for c in candidates:
        if c.exists():
            return c
    raise OSError("OmniChem shared library not found; build it first (scripts/build.sh or scripts\\build.ps1) "
                  "or set OMNICHEM_LIB. Looked in:\n  " + "\n  ".join(str(c) for c in candidates))


def _load() -> ctypes.CDLL:
    path = library_path()
    if platform.system() == "Windows" and hasattr(os, "add_dll_directory"):
        os.add_dll_directory(str(path.parent))  # dependent DLLs next to omnichem.dll
    return ctypes.CDLL(str(path))


_lib = _load()
_c_char_pp = ctypes.POINTER(ctypes.c_char_p)


def _sig(name, restype, *argtypes):
    fn = getattr(_lib, name)
    fn.restype = restype
    fn.argtypes = list(argtypes)
    return fn


_oc_free_string = _sig("oc_free_string", None, ctypes.c_void_p)
_oc_version = _sig("oc_version", ctypes.c_char_p)
_oc_parse = _sig("oc_parse_smiles_json", ctypes.c_void_p, ctypes.c_char_p)
_oc_canon = _sig("oc_canonical_smiles_json", ctypes.c_void_p, ctypes.c_char_p, ctypes.c_uint, ctypes.c_uint64)
_oc_kinetics = _sig("oc_kinetics_json", ctypes.c_void_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_double,
                    ctypes.c_double, _c_char_pp, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.c_double,
                    ctypes.c_size_t, ctypes.c_double)
_oc_reaction = _sig("oc_check_reaction_json", ctypes.c_void_p, ctypes.c_char_p, ctypes.c_double)
_oc_equil = _sig("oc_equilibrate_json", ctypes.c_void_p, _c_char_pp, ctypes.POINTER(ctypes.c_double),
                 ctypes.c_size_t, ctypes.c_double, ctypes.c_double)

_sys_new = _sig("oc_system_new", ctypes.c_void_p)
_sys_free = _sig("oc_system_free", None, ctypes.c_void_p)
_sys_add = _sig("oc_system_add_smiles", ctypes.c_int32, ctypes.c_void_p, ctypes.c_char_p)
_sys_set_moles = _sig("oc_system_set_moles", ctypes.c_int, ctypes.c_void_p, ctypes.c_uint32, ctypes.c_double)
_sys_set_cond = _sig("oc_system_set_conditions", ctypes.c_int, ctypes.c_void_p, ctypes.c_double, ctypes.c_double)
_sys_equil = _sig("oc_system_equilibrate", ctypes.c_int, ctypes.c_void_p)
_sys_atoms = _sig("oc_system_atom_count", ctypes.c_size_t, ctypes.c_void_p)
_sys_bonds = _sig("oc_system_bond_count", ctypes.c_size_t, ctypes.c_void_p)
_sys_species = _sig("oc_system_species_count", ctypes.c_size_t, ctypes.c_void_p)
_sys_struct = _sig("oc_system_structure_json", ctypes.c_void_p, ctypes.c_void_p)
_sys_thermo = _sig("oc_system_thermo_json", ctypes.c_void_p, ctypes.c_void_p)
_sys_diag = _sig("oc_system_diagnostics_json", ctypes.c_void_p, ctypes.c_void_p)

_columns = {
    "atom_element": (ctypes.c_uint8, "atoms"),
    "atom_charge": (ctypes.c_int8, "atoms"),
    "atom_implicit_h": (ctypes.c_uint8, "atoms"),
    "atom_x": (ctypes.c_double, "atoms"),
    "atom_y": (ctypes.c_double, "atoms"),
    "atom_z": (ctypes.c_double, "atoms"),
    "bond_begin": (ctypes.c_uint32, "bonds"),
    "bond_end": (ctypes.c_uint32, "bonds"),
    "bond_order": (ctypes.c_uint8, "bonds"),
    "moles": (ctypes.c_double, "species"),
}
_column_fns = {k: _sig(f"oc_system_{k}", ctypes.POINTER(t), ctypes.c_void_p) for k, (t, _) in _columns.items()}


def _take_json(ptr: int | None) -> dict | list:
    if not ptr:
        raise MemoryError("OmniChem returned NULL")
    try:
        return json.loads(ctypes.string_at(ptr).decode("utf-8"))
    finally:
        _oc_free_string(ptr)


def version() -> str:
    return _oc_version().decode()


def parse(smiles: str) -> dict:
    """Parse one SMILES; returns {'ok', 'structure', 'diagnostics'}."""
    return _take_json(_oc_parse(smiles.encode()))


def canonical(smiles: str, n_random: int = 0, seed: int = 1) -> dict:
    """Canonical SMILES, CIP labels, and optionally `n_random` randomised equivalent SMILES.

    Raises OmniChemError if the SMILES is invalid.
    """
    r = _take_json(_oc_canon(smiles.encode(), n_random, seed))
    if not r["ok"]:
        raise OmniChemError(r["diagnostics"])
    return r


def kinetics(mechanism: str, feed: dict[str, float], temperature_K: float, pressure_Pa: float = 101325.0,
             t_end: float = 1e-3, model: str = "adiabatic-v", samples: int = 50, rtol: float = 1e-8) -> dict:
    """Integrate a batch reactor. `mechanism` is .ocm text or a path to an .ocm file.

    Returns {'species', 'samples': [{'t','T','P','n'}], 'ignition_time', 'stats', 'diagnostics'}.
    Raises OmniChemError on invalid mechanisms or setups.
    """
    if "\n" not in mechanism and os.path.exists(mechanism):
        mechanism = Path(mechanism).read_text()
    names = (ctypes.c_char_p * len(feed))(*[k.encode() for k in feed])
    moles = (ctypes.c_double * len(feed))(*[float(v) for v in feed.values()])
    r = _take_json(_oc_kinetics(mechanism.encode(), model.encode(), temperature_K, pressure_Pa, names, moles,
                                len(feed), t_end, samples, rtol))
    if not r["ok"]:
        raise OmniChemError(r["diagnostics"])
    return r


def check_reaction(equation: str, temperature_K: float = 298.15) -> dict:
    """Validate stoichiometry (and compute ΔrG°, K when balanced)."""
    return _take_json(_oc_reaction(equation.encode(), temperature_K))


def equilibrate(species: dict[str, float] | Sequence[tuple[str, float]], temperature_K: float,
                pressure_Pa: float = 1e5) -> dict:
    """Gibbs-minimise an ideal-gas mixture. `species` maps SMILES -> initial moles."""
    items = list(species.items()) if isinstance(species, dict) else list(species)
    smi = (ctypes.c_char_p * len(items))(*[s.encode() for s, _ in items])
    mol = (ctypes.c_double * len(items))(*[float(n) for _, n in items])
    return _take_json(_oc_equil(smi, mol, len(items), temperature_K, pressure_Pa))


class System:
    """Owning handle around a ChemicalSystemState with zero-copy column access."""

    def __init__(self) -> None:
        self._h = _sys_new()
        if not self._h:
            raise MemoryError("oc_system_new failed")

    def close(self) -> None:
        if self._h:
            _sys_free(self._h)
            self._h = None

    def __del__(self) -> None:
        self.close()

    def __enter__(self) -> "System":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def add(self, smiles: str, moles: float = 0.0) -> int:
        idx = _sys_add(self._h, smiles.encode())
        if idx < 0:
            raise OmniChemError(self.diagnostics())
        if moles:
            _sys_set_moles(self._h, idx, moles)
        return idx

    def set_moles(self, species: int, moles: float) -> None:
        if not _sys_set_moles(self._h, species, moles):
            raise ValueError("invalid species index or amount")

    def set_conditions(self, temperature_K: float, pressure_Pa: float = 1e5) -> None:
        if not _sys_set_cond(self._h, temperature_K, pressure_Pa):
            raise ValueError("temperature and pressure must be positive")

    def equilibrate(self) -> bool:
        ok = bool(_sys_equil(self._h))
        if not ok:
            raise OmniChemError(self.diagnostics())
        return ok

    def diagnostics(self) -> list[dict]:
        return _take_json(_sys_diag(self._h))  # type: ignore[return-value]

    def structure(self) -> dict:
        return _take_json(_sys_struct(self._h))  # type: ignore[return-value]

    def thermo(self) -> dict:
        return _take_json(_sys_thermo(self._h))  # type: ignore[return-value]

    def _count(self, kind: str) -> int:
        return {"atoms": _sys_atoms, "bonds": _sys_bonds, "species": _sys_species}[kind](self._h)

    def column(self, name: str):
        """Zero-copy view of a SoA column (numpy array if numpy is installed, else ctypes array).

        The view is invalidated by the next add() on this system.
        """
        ctype, kind = _columns[name]
        n = self._count(kind)
        ptr = _column_fns[name](self._h)
        try:
            import numpy as np  # noqa: PLC0415
            return np.ctypeslib.as_array(ptr, shape=(n,)) if n else np.empty(0, dtype=np.dtype(ctype))
        except ImportError:
            return ctypes.cast(ptr, ctypes.POINTER(ctype * n)).contents if n else (ctype * 0)()

    def __getattr__(self, name: str):
        if name in _columns:
            return lambda: self.column(name)
        raise AttributeError(name)
