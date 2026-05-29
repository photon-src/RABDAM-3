"""Write simple Bnet reference-database CSV output.

Rejected entries are optional and written to a separate simple CSV for build
diagnostics. The accepted CSV is the file consumed by ``bnet.reference`` and
``bnet.percentile``.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import csv
from dataclasses import dataclass
import json
from pathlib import Path

from .process import (
    AcceptedBnetReferenceRow,
    PdbRedoProcessResult,
    RejectedBnetReferenceRow,
)


ACCEPTED_FIELDNAMES = (
    "",
    "PDB code",
    "Resolution (A)",
    "Bnet",
)

ACCEPTED_DETAILS_FIELDNAMES = (
    "PDB code",
    "Resolution (A)",
    "Bnet",
    "Rwork",
    "Rfree",
    "temperature_K",
    "wilson_b",
    "b_factor_restraint_weight",
    "bnet_site_count",
    "Asp/Glu carboxyl oxygen count",
    "asp_glu_residue_count",
    "median_bdamage",
    "left_area",
    "right_area",
    "atom_count",
    "non_hydrogen_atom_count",
    "protein_atom_count",
    "selected_atom_count",
    "bdamage_window_size",
    "has_protein",
    "has_nucleic_acid",
    "is_xray",
    "has_nonflat_protein_b_factors",
    "has_asp_glu_residue_with_total_occupancy_below_one",
    "experimental_methods",
    "metadata_warnings",
    "structure_check_warnings",
    "resolution_source",
    "r_work_source",
    "r_free_source",
    "temperature_source",
    "wilson_b_source",
    "b_factor_restraint_weight_source",
    "final_cif_path",
    "data_json_path",
)

REJECTED_FIELDNAMES = (
    "PDB code",
    "stage",
    "reason",
    "message",
    "Resolution (A)",
    "Rfree",
    "temperature_K",
    "Asp/Glu carboxyl oxygen count",
    "Bnet",
    "final_cif_path",
    "data_json_path",
    "exception_type",
    "traceback_text",
    "metadata_warnings",
    "structure_check_warnings",
)


@dataclass(frozen=True, slots=True)
class BnetDatabaseOutputPaths:
    """Output paths for a simple Bnet database build."""

    accepted_csv_path: Path
    accepted_details_csv_path: Path | None = None
    rejected_csv_path: Path | None = None


class BnetDatabaseOutputError(ValueError):
    """Raised when Bnet database output cannot be written."""


class BnetDatabaseCsvWriter:
    """Append accepted and rejected process results to simple CSV files."""

    def __init__(
        self,
        accepted_csv_path: str | Path,
        *,
        accepted_details_csv_path: str | Path | None = None,
        rejected_csv_path: str | Path | None = None,
        overwrite: bool = False,
    ) -> None:
        self.accepted_csv_path = Path(accepted_csv_path).expanduser()
        self.accepted_details_csv_path = (
            Path(accepted_details_csv_path).expanduser()
            if accepted_details_csv_path is not None
            else None
        )
        self.rejected_csv_path = (
            Path(rejected_csv_path).expanduser()
            if rejected_csv_path is not None
            else None
        )

        self.accepted_csv_path.parent.mkdir(parents=True, exist_ok=True)
        if self.accepted_details_csv_path is not None:
            self.accepted_details_csv_path.parent.mkdir(parents=True, exist_ok=True)
        if self.rejected_csv_path is not None:
            self.rejected_csv_path.parent.mkdir(parents=True, exist_ok=True)

        if overwrite:
            self.accepted_csv_path.unlink(missing_ok=True)
            if self.accepted_details_csv_path is not None:
                self.accepted_details_csv_path.unlink(missing_ok=True)
            if self.rejected_csv_path is not None:
                self.rejected_csv_path.unlink(missing_ok=True)

        _ensure_header(self.accepted_csv_path, ACCEPTED_FIELDNAMES)
        if self.accepted_details_csv_path is not None:
            _ensure_header(
                self.accepted_details_csv_path,
                ACCEPTED_DETAILS_FIELDNAMES,
            )
        if self.rejected_csv_path is not None:
            _ensure_header(self.rejected_csv_path, REJECTED_FIELDNAMES)

        self._accepted_count = _existing_data_row_count(self.accepted_csv_path)
        accepted_details_count = (
            _existing_data_row_count(self.accepted_details_csv_path)
            if self.accepted_details_csv_path is not None
            else self._accepted_count
        )
        if accepted_details_count != self._accepted_count:
            raise BnetDatabaseOutputError(
                "Accepted CSV and accepted-details CSV contain different "
                "numbers of data rows."
            )

        self._rejected_count = (
            _existing_data_row_count(self.rejected_csv_path)
            if self.rejected_csv_path is not None
            else 0
        )

    @property
    def accepted_count(self) -> int:
        """Return number of accepted rows written or already present."""

        return self._accepted_count

    @property
    def rejected_count(self) -> int:
        """Return number of rejected rows written or already present."""

        return self._rejected_count

    @property
    def paths(self) -> BnetDatabaseOutputPaths:
        """Return output paths."""

        return BnetDatabaseOutputPaths(
            accepted_csv_path=self.accepted_csv_path,
            accepted_details_csv_path=self.accepted_details_csv_path,
            rejected_csv_path=self.rejected_csv_path,
        )

    def write_result(self, result: PdbRedoProcessResult) -> None:
        """Write one accepted or rejected process result."""

        if result.accepted is not None:
            self.write_accepted(result.accepted)
            return

        if result.rejected is not None:
            self.write_rejected(result.rejected)
            return

        raise BnetDatabaseOutputError(
            "PdbRedoProcessResult contains neither accepted nor rejected row."
        )

    def write_accepted(self, row: AcceptedBnetReferenceRow) -> None:
        """Append one accepted Bnet database row in processing order."""

        output_row: dict[str, object] = {
            "": self._accepted_count,
            "PDB code": row.pdb_id.upper(),
            "Resolution (A)": _format_float(row.resolution_angstrom),
            "Bnet": _format_float(row.bnet),
        }

        _append_dict_row(self.accepted_csv_path, ACCEPTED_FIELDNAMES, output_row)
        if self.accepted_details_csv_path is not None:
            _append_dict_row(
                self.accepted_details_csv_path,
                ACCEPTED_DETAILS_FIELDNAMES,
                _accepted_details_dict(row),
            )
        self._accepted_count += 1

    def write_rejected(self, row: RejectedBnetReferenceRow) -> None:
        """Append one rejected build-log row, if rejected output is enabled."""

        if self.rejected_csv_path is None:
            return

        output_row = {
            "PDB code": row.pdb_id.upper(),
            "stage": row.stage.value,
            "reason": row.reason,
            "message": row.message,
            "Resolution (A)": _format_optional_float(row.resolution_angstrom),
            "Rfree": _format_optional_float(row.r_free),
            "temperature_K": _format_optional_float(row.temperature_k),
            "Asp/Glu carboxyl oxygen count": (
                "" if row.asp_glu_carboxyl_oxygen_count is None
                else str(row.asp_glu_carboxyl_oxygen_count)
            ),
            "Bnet": _format_optional_float(row.bnet),
            "final_cif_path": str(row.final_cif_path),
            "data_json_path": (
                "" if row.data_json_path is None else str(row.data_json_path)
            ),
            "exception_type": row.exception_type or "",
            "traceback_text": row.traceback_text or "",
            "metadata_warnings": _format_text_tuple(row.metadata_warnings),
            "structure_check_warnings": _format_text_tuple(
                row.structure_check_warnings
            ),
        }

        _append_dict_row(self.rejected_csv_path, REJECTED_FIELDNAMES, output_row)
        self._rejected_count += 1


def write_simple_bnet_database_csv(
    rows: Sequence[AcceptedBnetReferenceRow],
    path: str | Path,
    *,
    sort_by_bnet_descending: bool = True,
) -> None:
    """Write accepted rows to a simple Bnet reference CSV.

    This is useful for writing a completed database all at once. By default,
    rows are sorted by descending Bnet. For long builds, prefer
    ``BnetDatabaseCsvWriter`` so rows can be appended in processing order as
    each candidate finishes.
    """

    output_path = Path(path).expanduser()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    output_rows = list(rows)
    if sort_by_bnet_descending:
        output_rows.sort(key=lambda row: row.bnet, reverse=True)

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=ACCEPTED_FIELDNAMES)
        writer.writeheader()

        for index, row in enumerate(output_rows):
            writer.writerow(
                {
                    "": index,
                    "PDB code": row.pdb_id.upper(),
                    "Resolution (A)": _format_float(row.resolution_angstrom),
                    "Bnet": _format_float(row.bnet),
                }
            )


def _ensure_header(path: Path, fieldnames: tuple[str, ...]) -> None:
    if path.exists() and path.stat().st_size > 0:
        with path.open("r", newline="", encoding="utf-8") as handle:
            reader = csv.reader(handle)
            try:
                header = tuple(next(reader))
            except StopIteration:
                header = ()

        if header != fieldnames:
            raise BnetDatabaseOutputError(
                f"Existing CSV header does not match expected schema: {path}"
            )
        return

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()


def _accepted_details_dict(row: AcceptedBnetReferenceRow) -> dict[str, object]:
    return {
        "PDB code": row.pdb_id.upper(),
        "Resolution (A)": _format_float(row.resolution_angstrom),
        "Bnet": _format_float(row.bnet),
        "Rwork": _format_optional_float(row.r_work),
        "Rfree": _format_optional_float(row.r_free),
        "temperature_K": _format_optional_float(row.temperature_k),
        "wilson_b": _format_optional_float(row.wilson_b),
        "b_factor_restraint_weight": _format_optional_float(
            row.b_factor_restraint_weight
        ),
        "bnet_site_count": row.bnet_site_count,
        "Asp/Glu carboxyl oxygen count": row.asp_glu_carboxyl_oxygen_count,
        "asp_glu_residue_count": row.asp_glu_residue_count,
        "median_bdamage": _format_float(row.median_bdamage),
        "left_area": _format_float(row.left_area),
        "right_area": _format_float(row.right_area),
        "atom_count": row.atom_count,
        "non_hydrogen_atom_count": row.non_hydrogen_atom_count,
        "protein_atom_count": row.protein_atom_count,
        "selected_atom_count": row.selected_atom_count,
        "bdamage_window_size": row.bdamage_window_size,
        "has_protein": _format_bool(row.has_protein),
        "has_nucleic_acid": _format_bool(row.has_nucleic_acid),
        "is_xray": _format_bool(row.is_xray),
        "has_nonflat_protein_b_factors": _format_bool(
            row.has_nonflat_protein_b_factors
        ),
        "has_asp_glu_residue_with_total_occupancy_below_one": _format_bool(
            row.has_asp_glu_residue_with_total_occupancy_below_one
        ),
        "experimental_methods": _format_text_tuple(row.experimental_methods),
        "metadata_warnings": _format_text_tuple(row.metadata_warnings),
        "structure_check_warnings": _format_text_tuple(
            row.structure_check_warnings
        ),
        "resolution_source": row.resolution_source or "",
        "r_work_source": row.r_work_source or "",
        "r_free_source": row.r_free_source or "",
        "temperature_source": row.temperature_source or "",
        "wilson_b_source": row.wilson_b_source or "",
        "b_factor_restraint_weight_source": (
            row.b_factor_restraint_weight_source or ""
        ),
        "final_cif_path": "" if row.final_cif_path is None else str(row.final_cif_path),
        "data_json_path": "" if row.data_json_path is None else str(row.data_json_path),
    }


def _append_dict_row(
    path: Path,
    fieldnames: tuple[str, ...],
    row: Mapping[str, object],
) -> None:
    with path.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writerow(row)


def _existing_data_row_count(path: Path | None) -> int:
    if path is None or not path.exists() or path.stat().st_size == 0:
        return 0

    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        try:
            next(reader)
        except StopIteration:
            return 0

        return sum(1 for _ in reader)


def _format_float(value: float) -> str:
    return f"{value:.10g}"


def _format_optional_float(value: float | None) -> str:
    if value is None:
        return ""

    return _format_float(value)


def _format_bool(value: bool) -> str:
    return "true" if value else "false"


def _format_text_tuple(values: tuple[str, ...]) -> str:
    if not values:
        return ""

    return json.dumps(list(values), separators=(",", ":"))


__all__ = [
    "ACCEPTED_FIELDNAMES",
    "ACCEPTED_DETAILS_FIELDNAMES",
    "REJECTED_FIELDNAMES",
    "BnetDatabaseCsvWriter",
    "BnetDatabaseOutputError",
    "BnetDatabaseOutputPaths",
    "write_simple_bnet_database_csv",
]
