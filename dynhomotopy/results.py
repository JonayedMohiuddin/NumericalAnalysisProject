"""Result containers returned by the solvers and the hybrid pipeline."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class SolveResult:
    """Outcome of an iterative solver (NR, FDXB, GSH-NR).

    ``norms[0]`` is ||g(x)||_inf at the initial estimate and ``norms[i]`` the
    value after iteration i, which is how the paper tabulates the refinement.
    """

    method: str
    converged: bool
    iterations: int
    x: np.ndarray
    norms: list[float]
    factorizations: int = 0
    time: float = 0.0
    states: list[np.ndarray] | None = None  # x after every iteration (optional)


@dataclass
class Trajectory:
    """Discrete dynamic-homotopy pathway gamma(t): points (t_k, x_k)."""

    method: str
    K: float
    times: list[float]
    states: list[np.ndarray]
    norms: list[float]            # ||g(x_k)||_inf of the original PFP
    factorizations: int = 0
    time: float = 0.0
    failed: bool = False          # non-finite state or singular matrix

    @property
    def x_final(self) -> np.ndarray:
        """x(t)|_{t=1}: the initial estimate handed to NR / FDXB."""
        return self.states[-1]


@dataclass
class HybridResult:
    """Dynamic homotopy followed by a refining solver (NR or FDXB)."""

    trajectory: Trajectory
    refine: SolveResult | None
    time: float = 0.0
    extra: dict = field(default_factory=dict)

    @property
    def converged(self) -> bool:
        return self.refine is not None and self.refine.converged

    @property
    def factorizations(self) -> int:
        return self.trajectory.factorizations + (self.refine.factorizations if self.refine else 0)
