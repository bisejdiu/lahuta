from lahuta.pipeline import Pipeline


def test_pipeline_from_nmr_files_constructor() -> None:
    pipeline = Pipeline.from_nmr_files(["null.cif.gz"])
    assert isinstance(pipeline, Pipeline)
    # Graph inspection should not require the file to exist yet
    assert "Pipeline" in pipeline.describe()


def test_pipeline_from_md_trajectories_constructor() -> None:
    specs = [("null.gro", ["null.xtc"])]
    pipeline = Pipeline.from_md_trajectories(specs)
    assert isinstance(pipeline, Pipeline)
    assert "Pipeline" in pipeline.describe()
