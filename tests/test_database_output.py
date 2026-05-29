import csv
import json
from pathlib import Path
import tempfile
import unittest

from bnet.reference import load_bnet_reference_database
from database.output import (
    ACCEPTED_DETAILS_FIELDNAMES,
    ACCEPTED_FIELDNAMES,
    REJECTED_FIELDNAMES,
    BnetDatabaseCsvWriter,
    BnetDatabaseOutputError,
    write_simple_bnet_database_csv,
)
from database.process import (
    AcceptedBnetReferenceRow,
    PdbRedoProcessResult,
    PdbRedoProcessStage,
    RejectedBnetReferenceRow,
)


def make_accepted_row(
    pdb_id: str = "1abc",
    *,
    resolution_angstrom: float = 1.5,
    bnet: float = 1.2,
) -> AcceptedBnetReferenceRow:
    return AcceptedBnetReferenceRow(
        pdb_id=pdb_id,
        resolution_angstrom=resolution_angstrom,
        bnet=bnet,
        r_work=0.2,
        r_free=0.25,
        temperature_k=100.0,
        wilson_b=12.5,
        b_factor_restraint_weight=0.8,
        bnet_site_count=24,
        asp_glu_carboxyl_oxygen_count=24,
        asp_glu_residue_count=12,
        median_bdamage=1.1,
        left_area=0.4,
        right_area=0.6,
        atom_count=100,
        non_hydrogen_atom_count=90,
        protein_atom_count=80,
        selected_atom_count=70,
        bdamage_window_size=11,
        has_protein=True,
        has_nucleic_acid=False,
        is_xray=True,
        has_nonflat_protein_b_factors=True,
        has_asp_glu_residue_with_total_occupancy_below_one=False,
        experimental_methods=("X-RAY DIFFRACTION",),
        metadata_warnings=("metadata warning",),
        structure_check_warnings=("structure warning",),
        resolution_source="data_json:resolution",
        r_work_source="data_json:rwork",
        r_free_source="data_json:rfree",
        temperature_source="data_json:temperature",
        wilson_b_source="data_json:wilson_b",
        b_factor_restraint_weight_source="data_json:b_factor_restraint_weight",
        final_cif_path=Path("/tmp/1abc_final.cif"),
        data_json_path=Path("/tmp/data.json"),
    )


def make_rejected_row() -> RejectedBnetReferenceRow:
    return RejectedBnetReferenceRow(
        pdb_id="2def",
        stage=PdbRedoProcessStage.RABDAM,
        reason="rabdam_error",
        message="Calculation failed",
        final_cif_path=Path("/tmp/2def_final.cif"),
        data_json_path=Path("/tmp/data.json"),
        exception_type="ValueError",
        traceback_text="Traceback text",
        resolution_angstrom=2.0,
        r_free=0.3,
        temperature_k=100.0,
        asp_glu_carboxyl_oxygen_count=22,
        bnet=None,
        metadata_warnings=("metadata warning",),
        structure_check_warnings=("structure warning",),
    )


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


class BnetDatabaseOutputTests(unittest.TestCase):
    def test_writer_appends_minimal_details_and_rejected_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            accepted_path = root / "accepted.csv"
            details_path = root / "accepted_details.csv"
            rejected_path = root / "rejected.csv"

            writer = BnetDatabaseCsvWriter(
                accepted_path,
                accepted_details_csv_path=details_path,
                rejected_csv_path=rejected_path,
            )
            accepted_row = make_accepted_row()
            rejected_row = make_rejected_row()

            writer.write_result(
                PdbRedoProcessResult(
                    pdb_id=accepted_row.pdb_id,
                    accepted=accepted_row,
                )
            )
            writer.write_result(
                PdbRedoProcessResult(
                    pdb_id=rejected_row.pdb_id,
                    rejected=rejected_row,
                )
            )

            accepted_rows = read_csv_rows(accepted_path)
            detail_rows = read_csv_rows(details_path)
            rejected_rows = read_csv_rows(rejected_path)

        self.assertEqual(writer.accepted_count, 1)
        self.assertEqual(writer.rejected_count, 1)
        self.assertEqual(accepted_rows[0][""], "0")
        self.assertEqual(accepted_rows[0]["PDB code"], "1ABC")
        self.assertEqual(detail_rows[0]["PDB code"], "1ABC")
        self.assertEqual(detail_rows[0]["has_nonflat_protein_b_factors"], "true")
        self.assertEqual(
            json.loads(detail_rows[0]["metadata_warnings"]),
            ["metadata warning"],
        )
        self.assertEqual(rejected_rows[0]["PDB code"], "2DEF")
        self.assertEqual(rejected_rows[0]["traceback_text"], "Traceback text")
        self.assertEqual(
            json.loads(rejected_rows[0]["structure_check_warnings"]),
            ["structure warning"],
        )

    def test_accepted_output_round_trips_through_reference_loader(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            accepted_path = Path(temp_dir) / "accepted.csv"
            writer = BnetDatabaseCsvWriter(accepted_path)

            writer.write_accepted(make_accepted_row("1abc", bnet=1.2))
            writer.write_accepted(make_accepted_row("2def", bnet=2.4))

            database = load_bnet_reference_database(
                accepted_path,
                database_id="test_reference",
            )

        self.assertEqual(database.pdb_ids, ("1ABC", "2DEF"))
        self.assertEqual([entry.bnet for entry in database.entries], [1.2, 2.4])

    def test_writer_validates_existing_headers_before_appending(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            accepted_path = Path(temp_dir) / "accepted.csv"
            accepted_path.write_text(
                "pdb_id,resolution_angstrom,bnet\n"
                "1ABC,1.5,1.2\n",
                encoding="utf-8",
            )

            with self.assertRaises(BnetDatabaseOutputError):
                BnetDatabaseCsvWriter(accepted_path)

    def test_writer_rejects_mismatched_resume_counts(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            accepted_path = root / "accepted.csv"
            details_path = root / "accepted_details.csv"
            accepted_path.write_text(
                ",".join(ACCEPTED_FIELDNAMES)
                + "\n0,1ABC,1.5,1.2\n",
                encoding="utf-8",
            )
            details_path.write_text(
                ",".join(ACCEPTED_DETAILS_FIELDNAMES) + "\n",
                encoding="utf-8",
            )

            with self.assertRaises(BnetDatabaseOutputError):
                BnetDatabaseCsvWriter(
                    accepted_path,
                    accepted_details_csv_path=details_path,
                )

    def test_writer_overwrite_replaces_existing_files(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            accepted_path = root / "accepted.csv"
            rejected_path = root / "rejected.csv"
            accepted_path.write_text("bad\n", encoding="utf-8")
            rejected_path.write_text("bad\n", encoding="utf-8")

            writer = BnetDatabaseCsvWriter(
                accepted_path,
                rejected_csv_path=rejected_path,
                overwrite=True,
            )

            accepted_header = accepted_path.read_text(encoding="utf-8").splitlines()[0]
            rejected_header = rejected_path.read_text(encoding="utf-8").splitlines()[0]

        self.assertEqual(writer.accepted_count, 0)
        self.assertEqual(writer.rejected_count, 0)
        self.assertEqual(accepted_header, ",".join(ACCEPTED_FIELDNAMES))
        self.assertEqual(rejected_header, ",".join(REJECTED_FIELDNAMES))

    def test_rejected_output_can_be_disabled(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            accepted_path = Path(temp_dir) / "accepted.csv"
            writer = BnetDatabaseCsvWriter(accepted_path)

            writer.write_rejected(make_rejected_row())

        self.assertEqual(writer.rejected_count, 0)
        self.assertIsNone(writer.paths.rejected_csv_path)

    def test_batch_writer_sorts_by_bnet_by_default(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            accepted_path = Path(temp_dir) / "accepted.csv"

            write_simple_bnet_database_csv(
                (
                    make_accepted_row("1abc", bnet=1.2),
                    make_accepted_row("2def", bnet=2.4),
                ),
                accepted_path,
            )
            rows = read_csv_rows(accepted_path)
            database = load_bnet_reference_database(
                accepted_path,
                database_id="sorted_reference",
            )

        self.assertEqual([row[""] for row in rows], ["0", "1"])
        self.assertEqual(database.pdb_ids, ("2DEF", "1ABC"))
        self.assertEqual([entry.bnet for entry in database.entries], [2.4, 1.2])


if __name__ == "__main__":
    unittest.main()
