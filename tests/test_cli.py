from mesh_smoothing.cli import main


def test_main_reports_top_differences(capsys):
    rc = main(["--iterations", "1", "--top", "3", "--steps", "5"])
    out = capsys.readouterr().out
    assert rc == 0
    assert "Top 3 largest differences" in out


def test_main_saves_optimized_mesh(tmp_path):
    output = tmp_path / "out.vtk"
    rc = main(
        ["--iterations", "1", "--top", "1", "--steps", "5", "--save", str(output)]
    )
    assert rc == 0
    assert output.exists() and output.stat().st_size > 0
