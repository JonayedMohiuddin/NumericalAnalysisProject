import numpy as np
import scipy.sparse as sp
from pypower.dSbus_dV import dSbus_dV
from pypower.idx_bus import VA, VM
from pypower.idx_gen import GEN_BUS, VG

from ..problem import NonlinearProblem
from .case import bus_types
from .network import make_sbus, make_ybus


class PowerFlowProblem(NonlinearProblem):
    """Power balance equations (1)-(2) in polar form.

    The state is x = [theta_pv, theta_pq, V_pq] and the mismatch is
    g(x) = [dP_pv, dP_pq, dQ_pq] with dS = V * conj(Ybus V) - Sbus,
    the same ordering and sign as MATPOWER's newtonpf.
    """

    def __init__(self, case, ybus=None, sbus=None):
        self.case = case
        self.name = case.name
        if ybus is None:
            ybus = make_ybus(case.base_mva, case.bus, case.branch)
        self.ybus = sp.csr_matrix(ybus)
        self.sbus = make_sbus(case) if sbus is None else sbus
        self.ref, self.pv, self.pq = bus_types(case)
        self.pvpq = np.concatenate([self.pv, self.pq])
        self.npv, self.npq = len(self.pv), len(self.pq)

        # generator voltage set points hold |V| at PV and slack buses
        self.vm_fixed = np.ones(case.nb)
        self.vm_fixed[case.gen[:, GEN_BUS].astype(int)] = case.gen[:, VG]

    @property
    def n(self):
        return self.npv + 2 * self.npq

    def voltage(self, x):
        va = np.zeros(self.case.nb)
        vm = self.vm_fixed.copy()
        va[self.pvpq] = x[:self.npv + self.npq]
        vm[self.pq] = x[self.npv + self.npq:]
        return vm * np.exp(1j * va)

    def flat_start(self):
        return np.concatenate([np.zeros(self.npv + self.npq), np.ones(self.npq)])

    def case_start(self):
        """Initial guess stored in the case file, with the slack angle moved to 0."""
        va = np.deg2rad(self.case.bus[:, VA])
        va -= va[self.ref[0]]
        vm = self.case.bus[:, VM]
        return np.concatenate([va[self.pvpq], vm[self.pq]])

    def angle_index(self, bus):
        return int(np.flatnonzero(self.pvpq == bus)[0])

    def magnitude_index(self, bus):
        return self.npv + self.npq + int(np.flatnonzero(self.pq == bus)[0])

    def complex_mismatch(self, v):
        return v * np.conj(self.ybus @ v) - self.sbus

    def g(self, x):
        mis = self.complex_mismatch(self.voltage(x))
        return np.concatenate([mis[self.pvpq].real, mis[self.pq].imag])

    def jacobian(self, x):
        ds_dvm, ds_dva = dSbus_dV(self.ybus, self.voltage(x))
        ds_dvm, ds_dva = sp.csr_matrix(ds_dvm), sp.csr_matrix(ds_dva)

        j11 = ds_dva[self.pvpq][:, self.pvpq].real
        j12 = ds_dvm[self.pvpq][:, self.pq].real
        j21 = ds_dva[self.pq][:, self.pvpq].imag
        j22 = ds_dvm[self.pq][:, self.pq].imag
        return sp.bmat([[j11, j12], [j21, j22]], format="csc")
