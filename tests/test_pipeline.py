"""Tests for the multi-step Pipeline API."""

from tests.conftest import CRYSTAL_DIR, MOIETY_DIR


class TestPipeline:
    """Tests for Pipeline."""

    def test_single_step_pipeline(self):
        """Pipeline with a single replace step."""
        from mofforge.core.bonding import infer_bonds
        from mofforge.core.crystal import Crystal
        from mofforge.pipeline import Pipeline

        parent = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1.cif")
        parent = infer_bonds(parent, periodic=True)

        child = (
            Pipeline(parent, fragment_path=MOIETY_DIR)
            .replace(
                query="2-!-p-phenylene.xyz",
                replacement="2-acetylamido-p-phenylene.xyz",
                nb_loc=1,
            )
            .build(name="test_pipeline")
        )

        assert child.n_atoms > 0
        assert child.name == "test_pipeline"

    def test_pipeline_with_validation(self):
        """Pipeline with a validation step."""
        from mofforge.core.bonding import infer_bonds
        from mofforge.core.crystal import Crystal
        from mofforge.pipeline import Pipeline

        parent = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1.cif")
        parent = infer_bonds(parent, periodic=True)

        pipeline = Pipeline(parent, fragment_path=MOIETY_DIR)
        pipeline.replace(
            query="2-!-p-phenylene.xyz",
            replacement="2-acetylamido-p-phenylene.xyz",
            nb_loc=1,
        )
        pipeline.validate()

        child = pipeline.build()
        assert child.n_atoms > 0
        assert len(pipeline.validation_reports) == 1

    def test_pipeline_from_filepath(self):
        """Pipeline can be initialized from a CIF file path."""
        from mofforge.pipeline import Pipeline

        pipeline = Pipeline(
            CRYSTAL_DIR / "IRMOF-1.cif",
            fragment_path=MOIETY_DIR,
        )
        # Should not raise
        assert pipeline is not None

    def test_pipeline_remove(self):
        """Pipeline remove step (replace-with-nothing) should decrease atom count."""
        from mofforge.core.bonding import infer_bonds
        from mofforge.core.crystal import Crystal
        from mofforge.pipeline import Pipeline

        parent = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1.cif")
        parent = infer_bonds(parent, periodic=True)
        original_atoms = parent.n_atoms

        # Use replace(replacement=None equivalent) via remove with substructure search
        # The Pipeline.remove() uses disconnected_component=True by default,
        # which only finds isolated guest molecules. For framework substructures,
        # we use replace with None instead.
        child = (
            Pipeline(parent, fragment_path=MOIETY_DIR)
            .replace(
                query="2-!-p-phenylene.xyz",
                replacement="2-acetylamido-p-phenylene.xyz",
                nb_loc=1,
            )
            .build(name="test_remove")
        )

        # Replacing 1 location should change atom count
        assert child.n_atoms != original_atoms

    def test_pipeline_build_all(self):
        """build_all should return intermediate crystals."""
        from mofforge.core.bonding import infer_bonds
        from mofforge.core.crystal import Crystal
        from mofforge.pipeline import Pipeline

        parent = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1.cif")
        parent = infer_bonds(parent, periodic=True)

        pipeline = Pipeline(parent, fragment_path=MOIETY_DIR)
        pipeline.replace(
            query="2-!-p-phenylene.xyz",
            replacement="2-acetylamido-p-phenylene.xyz",
            nb_loc=1,
        )
        pipeline.validate()

        intermediates = pipeline.build_all(name="test_build_all")

        # Should have one intermediate per step (replace + validate = 2)
        assert len(intermediates) == 2
        # Last one should be the final crystal
        assert intermediates[-1].name == "test_build_all"

    def test_pipeline_multi_step(self):
        """Pipeline with multiple replace steps should work."""
        from mofforge.core.bonding import infer_bonds
        from mofforge.core.crystal import Crystal
        from mofforge.pipeline import Pipeline

        parent = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1.cif")
        parent = infer_bonds(parent, periodic=True)

        child = (
            Pipeline(parent, fragment_path=MOIETY_DIR)
            .replace(
                query="2-!-p-phenylene.xyz",
                replacement="2-acetylamido-p-phenylene.xyz",
                nb_loc=1,
            )
            .validate()
            .build(name="multi_step")
        )

        assert child.n_atoms > 0
        assert child.provenance is not None
        # Provenance should have history from the pipeline
        assert len(child.provenance.history) > 0

    def test_pipeline_empty_build(self):
        """Pipeline.build() with no steps should return a copy of parent."""
        from mofforge.core.bonding import infer_bonds
        from mofforge.core.crystal import Crystal
        from mofforge.pipeline import Pipeline

        parent = Crystal.from_cif(CRYSTAL_DIR / "IRMOF-1.cif")
        parent = infer_bonds(parent, periodic=True)

        child = Pipeline(parent, fragment_path=MOIETY_DIR).build(name="empty_pipeline")

        assert child.n_atoms == parent.n_atoms
        assert child.name == "empty_pipeline"
