from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy.io as sio
from pypower.bustypes import bustypes
from pypower.ext2int import ext2int


@dataclass
class Case:
    """A MATPOWER case with buses renumbered 0..nb-1 (internal indexing)."""

    name: str
    base_mva: float
    bus: np.ndarray
    gen: np.ndarray
    branch: np.ndarray
    ext_bus_ids: np.ndarray

    @property
    def nb(self):
        return self.bus.shape[0]

    def internal_bus(self, ext_id):
        hits = np.flatnonzero(self.ext_bus_ids == ext_id)
        if hits.size == 0:
            raise KeyError(f"bus {ext_id} not in case {self.name}")
        return int(hits[0])


def load_case(path):
    """Read a MATPOWER .mat case and convert it to internal indexing with PYPOWER's ext2int."""
    path = Path(path)
    data = sio.loadmat(path)
    mpc = data["mpc"] if "mpc" in data else data["mpc_m"]
    mpc = mpc[0, 0]
    ppc = {
        "baseMVA": float(np.asarray(mpc["baseMVA"], dtype=float).ravel()[0]),
        "bus": np.asarray(mpc["bus"], dtype=float),
        "gen": np.asarray(mpc["gen"], dtype=float),
        "branch": np.asarray(mpc["branch"], dtype=float),
    }
    ppc = ext2int(ppc)
    ext_ids = ppc["order"]["bus"]["i2e"].astype(int)
    return Case(path.stem, ppc["baseMVA"], ppc["bus"], ppc["gen"], ppc["branch"], ext_ids)


def bus_types(case):
    """Indices of the reference, PV and PQ buses."""
    return bustypes(case.bus, case.gen)
