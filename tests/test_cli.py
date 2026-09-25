from extensions.__main__ import main as extensions_main
from improvements.__main__ import main as improvements_main


def test_improvements_command(capsys):
    improvements_main(["case69limit", "--corrector", "--multiplier"])
    out = capsys.readouterr().out
    assert "paper" in out and "OM + PC" in out and "solved" in out


def test_extensions_commands(capsys):
    extensions_main(["spectrum", "case69limit"])
    extensions_main(["scratch", "case69limit", "--refiner", "FDXB"])
    out = capsys.readouterr().out
    assert "inverse iteration" in out and "our LU 8, SuperLU 8" in out
