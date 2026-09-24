"""A wider test bed: the paper's 9 cases plus 25 more from the same two Zenodo records."""

from dynhomotopy import datasets
from dynhomotopy.powerflow import PowerFlowProblem, load_case

EXTRA_ILL = ["case6024", "case6243", "case6748", "case7092", "case9961", "case10595", "case12110"]
EXTRA_LIMIT = [
    "case14limit", "case_ieee30limit", "case57limit", "case89pegaselimit", "case118limit",
    "case300limit", "case1354pegaselimit", "case2383wplimit", "case2736splimit",
    "case2737soplimit", "case2746woplimit", "case2746wplimit", "case2869pegaselimit",
    "case3012wplimit", "case3120splimit", "case3375wplimit", "case9241pegaselimit",
    "case13659pegaselimit",
]
EXTRA = EXTRA_ILL + EXTRA_LIMIT
ALL = datasets.ALL + EXTRA

_loaded = {}


def load(name):
    if name in datasets.SYSTEMS:
        return datasets.load_problem(name)
    if name not in _loaded:
        record = datasets.ILL_CONDITIONED_RECORD if name in EXTRA_ILL else datasets.LIMIT_CASES_RECORD
        _loaded[name] = PowerFlowProblem(load_case(datasets.download(name, record)))
    return _loaded[name]
