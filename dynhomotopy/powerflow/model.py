"""The power-flow problem (PFP) as a NonlinearProblem (Section 2, eqs. 1-3).

State vector (MATPOWER ordering):  x = [theta_pv; theta_pq; V_pq]  (radians, p.u.)
Mismatch:                          g(x) = [dP_pv; dP_pq; dQ_pq],
                                   dS = V .* conj(Ybus V) - Sbus

The sign convention is MATPOWER's (calculated minus specified), under which the
Jacobian has a positive-dominant diagonal, so the "easy" homotopy term
(1 - t) K I of eq. (23) acts as a regularising shift.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from ..problem import NonlinearProblem
from . import idx
from .case import Case, bus_types
from .network import make_sbus, make_ybus


class PowerFlowProblem(NonlinearProblem):
    """Polar power-balance equations of a MATPOWER case.

    ``ybus`` and ``sbus`` may be overridden, which the GSH-NR static homotopy
    uses to solve modified networks with the same bus classification.
    """

    def __init__(self, case: Case, ybus: sp.spmatrix | None = None,
                 sbus: np.ndarray | None = None):
        self.case = case
        self.name = case.name
        self.ybus = sp.csr_matrix(make_ybus(case.base_mva, case.bus, case.branch) if ybus is None else ybus)
        self.sbus = make_sbus(case) if sbus is None else sbus
        self.ref, self.pv, self.pq = bus_types(case)
        self.pvpq = np.concatenate([self.pv, self.pq])
        self.npv, self.npq = self.pv.size, self.pq.size

        # Voltage magnitudes held by generators (PV and slack buses), as in MATPOWER.
        vm = np.ones(case.nb)
        gbus = case.gen[:, idx.GEN_BUS].astype(int)
        vm[gbus] = case.gen[:, idx.VG]
        self._vm_fixed = vm

    # ---- state <-> voltage mapping ---------------------------------------------
    @property
    def n(self) -> int:
        return self.npv + 2 * self.npq

    def voltage(self, x: np.ndarray) -> np.ndarray:
        """Complex bus voltages for state x (slack angle is the 0 reference)."""
        va = np.zeros(self.case.nb)
        vm = self._vm_fixed.copy()
        va[self.pvpq] = x[: self.npv + self.npq]
        vm[self.pq] = x[self.npv + self.npq:]
        return vm * np.exp(1j * va)

    def state(self, v: np.ndarray) -> np.ndarray:
        """State vector of a complex voltage profile."""
        return np.concatenate([np.angle(v[self.pvpq]), np.abs(v[self.pq])])

    def flat_start(self) -> np.ndarray:
        """Flat start: all angles 0, PQ magnitudes 1 p.u. (Section 3.1)."""
        return np.concatenate([np.zeros(self.npv + self.npq), np.ones(self.npq)])

    def case_start(self) -> np.ndarray:
        """Native MATPOWER initial estimate (the reference NR-MAT run).

        Angles are shifted so the slack angle is 0, which leaves the problem
        unchanged (only angle differences enter g).
        """
        va = np.deg2rad(self.case.bus[:, idx.VA])
        va = va - va[self.ref[0]]
        vm = self.case.bus[:, idx.VM]
        return np.concatenate([va[self.pvpq], vm[self.pq]])

    def angle_index(self, bus: int) -> int:
        """Position of a bus angle in x."""
        return int(np.flatnonzero(self.pvpq == bus)[0])

    def magnitude_index(self, bus: int) -> int:
        """Position of a PQ-bus magnitude in x."""
        return self.npv + self.npq + int(np.flatnonzero(self.pq == bus)[0])

    # ---- equations ----------------------------------------------------------------
    def complex_mismatch(self, v: np.ndarray) -> np.ndarray:
        return v * np.conj(self.ybus @ v) - self.sbus

    def g(self, x: np.ndarray) -> np.ndarray:
        mis = self.complex_mismatch(self.voltage(x))
        return np.concatenate([mis[self.pvpq].real, mis[self.pq].imag])

    def jacobian(self, x: np.ndarray) -> sp.csc_matrix:
        """Polar Jacobian from MATPOWER's ``dSbus_dV``."""
        v = self.voltage(x)
        y = self.ybus
        ibus = y @ v
        diag_v = sp.diags(v)
        diag_vnorm = sp.diags(v / np.abs(v))
        ds_dvm = diag_v @ np.conj(y @ diag_vnorm) + sp.diags(np.conj(ibus)) @ diag_vnorm
        ds_dva = 1j * diag_v @ np.conj(sp.diags(ibus) - y @ diag_v)

        ds_dva = sp.csr_matrix(ds_dva)
        ds_dvm = sp.csr_matrix(ds_dvm)
        rows_p, rows_q = ds_dva[self.pvpq], ds_dva[self.pq]
        j11 = rows_p[:, self.pvpq].real
        j21 = rows_q[:, self.pvpq].imag
        rows_p, rows_q = ds_dvm[self.pvpq], ds_dvm[self.pq]
        j12 = rows_p[:, self.pq].real
        j22 = rows_q[:, self.pq].imag
        return sp.bmat([[j11, j12], [j21, j22]], format="csc")
