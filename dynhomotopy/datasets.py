"""Registry of the paper's test systems (Section 4.1, Table 1) and their download.

Ill-conditioned cases: Zenodo record 3514739, "Matpower Ill-conditioned systems" [46]
Stressed (limit) cases: Zenodo record 3491654, "Matpower Limit cases" [47]
"""

from __future__ import annotations

import urllib.request
from dataclasses import dataclass
from pathlib import Path

from .powerflow import PowerFlowProblem, load_case

DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "cases"


@dataclass(frozen=True)
class TestSystem:
    name: str             # file stem
    label: str            # name printed in the paper's tables
    short: str            # mnemonic used in the paper's figure legends
    kind: str             # "ill" (ill-conditioned) or "limit" (stressed)
    zenodo: str           # Zenodo record id
    paper_base_mva: float | None = None  # see note below


# case69limit and case141limit are stored with baseMVA = 10. The mismatch
# norms the paper reports for them in Tables 2-5 (0.012 and 0.006 at the flat
# start) and their BE(NR) = 1 iteration in Table 7 are reproduced only on a
# 100 MVA base; the other cases already use 100 MVA.
SYSTEMS: dict[str, TestSystem] = {s.name: s for s in [
    TestSystem("case18482", "case18482", "18k", "ill", "3514739"),
    TestSystem("case27318", "case27318", "27k", "ill", "3514739"),
    TestSystem("case36964", "case36964", "36k", "ill", "3514739"),
    TestSystem("case54636", "case54636", "54k", "ill", "3514739"),
    TestSystem("case109272", "case109272", "109k", "ill", "3514739"),
    TestSystem("case69limit", "case69limit", "69", "limit", "3491654", paper_base_mva=100.0),
    TestSystem("case141limit", "case141limit", "141", "limit", "3491654", paper_base_mva=100.0),
    TestSystem("case_ACTIVSg500limit", "case500limit", "500", "limit", "3491654"),
    TestSystem("case_ACTIVSg2000limit", "case2000limit", "lim", "limit", "3491654"),
]}

ILL_CONDITIONED = [s for s, v in SYSTEMS.items() if v.kind == "ill"]
LIMIT = [s for s, v in SYSTEMS.items() if v.kind == "limit"]
ALL = ILL_CONDITIONED + LIMIT


def case_path(name: str) -> Path:
    """Local path of a case, downloading it from Zenodo on first use."""
    path = DATA_DIR / f"{name}.mat"
    if not path.exists():
        system = SYSTEMS[name]
        url = f"https://zenodo.org/api/records/{system.zenodo}/files/{name}.mat/content"
        DATA_DIR.mkdir(parents=True, exist_ok=True)
        print(f"downloading {name} from {url}")
        urllib.request.urlretrieve(url, path)
    return path


_cache: dict[tuple[str, bool], PowerFlowProblem] = {}


def load_problem(name: str, paper_base: bool = True) -> PowerFlowProblem:
    """PowerFlowProblem of a test system (cached).

    ``paper_base=False`` keeps the base MVA stored in the file.
    """
    key = (name, paper_base)
    if key not in _cache:
        case = load_case(case_path(name))
        base = SYSTEMS[name].paper_base_mva if name in SYSTEMS else None
        if paper_base and base is not None:
            case.base_mva = base
        _cache[key] = PowerFlowProblem(case)
    return _cache[key]
