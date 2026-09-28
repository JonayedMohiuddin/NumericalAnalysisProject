from improvements.__main__ import main


def test_improvements_command(capsys):
    main(["case69limit", "--corrector", "--multiplier"])
    out = capsys.readouterr().out
    assert "paper" in out and "OM + PC" in out and "solved" in out


def test_study_commands(capsys):
    main(["spectrum", "case69limit"])
    main(["scratch", "case69limit", "--refiner", "FDXB"])
    out = capsys.readouterr().out
    assert "inverse iteration" in out and "our LU 8, SuperLU 8" in out
