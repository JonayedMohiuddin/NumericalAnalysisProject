from dataclasses import dataclass

import numpy as np


@dataclass
class SolveResult:
    # Result of NR, FDXB or GSH-NR. norms[0] is the mismatch at the initial guess

    method: str
    converged: bool
    iterations: int
    x: np.ndarray
    norms: list
    factorizations: int = 0
    time: float = 0.0
    states: list | None = None


@dataclass
class Trajectory:
    # Points (t_k, x_k) computed along the homotopy path

    method: str
    K: float
    times: list
    states: list
    norms: list
    factorizations: int = 0
    time: float = 0.0
    failed: bool = False

    @property
    def x_final(self):
        return self.states[-1]


@dataclass
class HybridResult:
    trajectory: Trajectory
    refine: SolveResult | None
    time: float = 0.0

    @property
    def converged(self):
        return self.refine is not None and self.refine.converged

    @property
    def factorizations(self):
        refine_lus = self.refine.factorizations if self.refine else 0
        return self.trajectory.factorizations + refine_lus
