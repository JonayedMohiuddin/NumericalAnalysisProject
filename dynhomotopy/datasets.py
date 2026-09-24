"""The test systems of Section 4.1, downloaded from Zenodo when first needed."""

import urllib.request
from dataclasses import dataclass
from pathlib import Path

from .powerflow import PowerFlowProblem, load_case

DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "cases"

ILL_CONDITIONED_RECORD = "3514739"
LIMIT_CASES_RECORD = "3491654"


@dataclass(frozen=True)
class System:
    name: str
    label: str
    short: str
    kind: str
    zenodo: str
    base_mva: float | None = None


# case69limit and case141limit are stored with baseMVA = 10, but the paper's
# numbers for them in Tables 2-5 only come out on a 100 MVA base.
SYSTEMS = {s.name: s for s in [
    System("case18482", "case18482", "18k", "ill", ILL_CONDITIONED_RECORD),
    System("case27318", "case27318", "27k", "ill", ILL_CONDITIONED_RECORD),
    System("case36964", "case36964", "36k", "ill", ILL_CONDITIONED_RECORD),
    System("case54636", "case54636", "54k", "ill", ILL_CONDITIONED_RECORD),
    System("case109272", "case109272", "109k", "ill", ILL_CONDITIONED_RECORD),
    System("case69limit", "case69limit", "69", "limit", LIMIT_CASES_RECORD, base_mva=100),
    System("case141limit", "case141limit", "141", "limit", LIMIT_CASES_RECORD, base_mva=100),
    System("case_ACTIVSg500limit", "case500limit", "500", "limit", LIMIT_CASES_RECORD),
    System("case_ACTIVSg2000limit", "case2000limit", "lim", "limit", LIMIT_CASES_RECORD),
]}

ILL_CONDITIONED = [s.name for s in SYSTEMS.values() if s.kind == "ill"]
LIMIT = [s.name for s in SYSTEMS.values() if s.kind == "limit"]
ALL = ILL_CONDITIONED + LIMIT


def case_path(name):
    path = DATA_DIR / f"{name}.mat"
    if not path.exists():
        url = f"https://zenodo.org/api/records/{SYSTEMS[name].zenodo}/files/{name}.mat/content"
        DATA_DIR.mkdir(parents=True, exist_ok=True)
        print(f"downloading {name} from {url}")
        urllib.request.urlretrieve(url, path)
    return path


_loaded = {}


def load_problem(name, paper_base=True):
    """Load a test system. paper_base=False keeps the base MVA stored in the file."""
    key = (name, paper_base)
    if key not in _loaded:
        case = load_case(case_path(name))
        if paper_base and SYSTEMS[name].base_mva:
            case.base_mva = SYSTEMS[name].base_mva
        _loaded[key] = PowerFlowProblem(case)
    return _loaded[key]
