from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy.io as sio

from . import idx


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
    path = Path(path)
    data = sio.loadmat(path)
    mpc = data["mpc"] if "mpc" in data else data["mpc_m"]
    mpc = mpc[0, 0]
    base_mva = float(np.asarray(mpc["baseMVA"], dtype=float).ravel()[0])
    bus = np.asarray(mpc["bus"], dtype=float)
    gen = np.asarray(mpc["gen"], dtype=float)
    branch = np.asarray(mpc["branch"], dtype=float)
    return ext2int(path.stem, base_mva, bus, gen, branch)


def ext2int(name, base_mva, bus, gen, branch):
    """Same as MATPOWER's ext2int: drop isolated buses and offline elements, renumber buses."""
    ext_ids = bus[:, idx.BUS_I].astype(int)
    row_of = {b: i for i, b in enumerate(ext_ids)}
    bus_on = bus[:, idx.BUS_TYPE] != idx.NONE

    def connected(bus_numbers):
        return np.array([bus_on[row_of[int(b)]] for b in bus_numbers], dtype=bool)

    gen_on = (gen[:, idx.GEN_STATUS] > 0) & connected(gen[:, idx.GEN_BUS])
    br_on = ((branch[:, idx.BR_STATUS] > 0) & connected(branch[:, idx.F_BUS])
             & connected(branch[:, idx.T_BUS]))

    bus = bus[bus_on].copy()
    gen = gen[gen_on].copy()
    branch = branch[br_on].copy()
    ext_ids = ext_ids[bus_on]
    new_index = {b: i for i, b in enumerate(ext_ids)}

    bus[:, idx.BUS_I] = np.arange(len(bus))
    for matrix, cols in ((gen, [idx.GEN_BUS]), (branch, [idx.F_BUS, idx.T_BUS])):
        for c in cols:
            matrix[:, c] = [new_index[int(b)] for b in matrix[:, c]]

    # MATPOWER sorts generators by bus; this decides which VG is used
    # when a bus has more than one generator.
    gen = gen[np.argsort(gen[:, idx.GEN_BUS], kind="stable")]
    return Case(name, base_mva, bus, gen, branch, ext_ids)


def bus_types(case):
    """Indices of the reference, PV and PQ buses (MATPOWER's bustypes)."""
    has_gen = np.zeros(case.nb, dtype=bool)
    has_gen[case.gen[:, idx.GEN_BUS].astype(int)] = True
    btype = case.bus[:, idx.BUS_TYPE]
    ref = np.flatnonzero((btype == idx.REF) & has_gen)
    pv = np.flatnonzero((btype == idx.PV) & has_gen)
    pq = np.flatnonzero((btype == idx.PQ) | ~has_gen)
    if ref.size == 0:
        ref, pv = pv[:1], pv[1:]
    return ref, pv, pq
