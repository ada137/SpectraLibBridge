
import argparse
import json
from pathlib import Path
from datetime import date, datetime
import math
import tempfile


def validated_peaks(record):
    """Return sorted peak pairs; reject empty, malformed, or nonfinite data."""
    raw_peaks = record["peaks"]
    if not isinstance(raw_peaks, list) or not raw_peaks:
        raise ValueError("Spectrum must contain a nonempty list of peaks")
    peaks = []
    for pair in raw_peaks:
        if not isinstance(pair, (list, tuple)) or len(pair) != 2:
            raise ValueError("Each peak must contain exactly m/z and intensity")
        if any(isinstance(value, bool) for value in pair):
            raise ValueError("Peak values cannot be booleans")
        mz, intensity = map(float, pair)
        if not math.isfinite(mz) or not math.isfinite(intensity) or mz <= 0 or intensity < 0:
            raise ValueError("Peaks must have positive finite m/z and nonnegative finite intensity")
        peaks.append((mz, intensity))
    return sorted(peaks)


def load_json_records(path_to_json_file):
    """Read MZmine JSON Lines, yielding (source line number, record) pairs."""
    with Path(path_to_json_file).open(encoding="utf-8-sig") as source:
        for line_number, line in enumerate(source, 1):
            if not line.strip():
                continue
            try:
                record = json.loads(line)
                if not isinstance(record, dict):
                    raise ValueError("Expected a JSON object for each spectrum")
                validated_peaks(record)
            except (ValueError, KeyError, TypeError, OverflowError) as error:
                raise ValueError(f"Input line {line_number}: {error}") from error
            yield line_number, record

def text(value):
    if isinstance(value, (list, dict)):
        value = json.dumps(value, ensure_ascii=False)
    return str(value).replace("\r", " ").replace("\n", " ")


def convert_record_to_txt(record, number):
    """
    Convert one source spectrum into MassBank-style TXT content.

    Sort and validate the peaks, map source metadata to TXT fields,
    and calculate relative intensities on a scale from 0 to 999.
    Preserve the original metadata in a SOURCE_METADATA comment.

    Parameters
    ----------
    record : dict
        One spectrum loaded from the source JSON.
        Must contain "peaks" as a list of [m/z, intensity] pairs.
        Other fields, such as "compound_name", "formula", and
        "collision_energy", are optional.
    number : int
        Spectrum number used to generate the local accession.
        For example, 1 produces "LOCAL000001".

    Returns
    -------
    accession : str
        Generated identifier, also used for the output filename.
    content : str
        Complete TXT record, including metadata, the peak table,
        and the closing "//" marker.

    Raises
    ------
    KeyError
        If the required "peaks" field is missing.
    ValueError
        If peaks are empty, contain invalid numbers, have m/z
        values less than or equal to zero, or negative intensities.

    Notes
    -----
    Does not modify the source record or write a file.
    Uses the original source data, not a matchms Spectrum object.
    The output has not been validated for MassBank submission.
    """
    peaks = validated_peaks(record)
    maximum = max(i for _, i in peaks)
    name = record.get("compound_name") or "Unknown"
    polarity = record.get("polarity")
    ion_mode = {"+": "POSITIVE", "-": "NEGATIVE"}.get(polarity)
    accession = f"LOCAL{number:06d}"
    title = [text(name)]
    if record.get("instrument"):
        title.append(text(record["instrument"]))
    if record.get("ms_level") is not None:
        title.append(f"MS{record['ms_level']}")
    if record.get("collision_energy") is not None:
        title.append(f"CE: {text(record['collision_energy'])} (units unspecified)")
    if record.get("adduct"):
        title.append(text(record["adduct"]))
    lines = [
        f"ACCESSION: {accession}",
        f"RECORD_TITLE: {'; '.join(title)}",
        f"DATE: {date.today():%Y.%m.%d}",
        "COMMENT: Local MassBank-style export; accession and date generated during conversion.",
        "COMMENT: Source licence not supplied; record has not been validated for MassBank submission.",
    ]

    def add(label, value):
        if value is not None and value != "":
            lines.append(f"{label}: {text(value)}")

    # Preserve source metadata, including fields without a direct target mapping.
    metadata = {key: value for key, value in record.items() if key != "peaks"}
    add("COMMENT", "SOURCE_METADATA " + json.dumps(metadata, ensure_ascii=False))
    add("CH$NAME", name)
    add("CH$FORMULA", record.get("formula"))
    add("CH$EXACT_MASS", record.get("exact_mass"))
    add("CH$SMILES", record.get("smiles"))
    # The reference's CH$IUPAC field contains an InChI, not an IUPAC name.
    add("CH$IUPAC", record.get("inchi"))
    if record.get("inchikey"):
        add("CH$LINK", "INCHIKEY " + text(record["inchikey"]))
    add("AC$INSTRUMENT", record.get("instrument"))
    if record.get("ms_level") is not None:
        add("AC$MASS_SPECTROMETRY", f"MS_TYPE MS{record['ms_level']}")
    if ion_mode:
        add("AC$MASS_SPECTROMETRY", "ION_MODE " + ion_mode)
    if record.get("ion_source"):
        add("AC$MASS_SPECTROMETRY", "IONIZATION " + text(record["ion_source"]))
    fragmentation = record.get("fragmentation_method") or record.get("fragmenation_method")
    if fragmentation:
        add("AC$MASS_SPECTROMETRY", "FRAGMENTATION_MODE " + text(fragmentation))
    if record.get("collision_energy") is not None:
        add("AC$MASS_SPECTROMETRY", "COLLISION_ENERGY "
            + text(record["collision_energy"]) + " (units unspecified in source)")
    if record.get("rt") is not None:
        add("AC$CHROMATOGRAPHY", "RETENTION_TIME " + text(record["rt"])
            + " (units unspecified in source)")
    if record.get("precursor_mz") is not None:
        add("MS$FOCUSED_ION", "PRECURSOR_M/Z " + text(record["precursor_mz"]))
    if record.get("adduct"):
        add("MS$FOCUSED_ION", "PRECURSOR_TYPE " + text(record["adduct"]))
    if record.get("merge_type"):
        add("MS$DATA_PROCESSING", "MERGE_TYPE " + text(record["merge_type"]))
    lines.extend([f"PK$NUM_PEAK: {len(peaks)}", "PK$PEAK: m/z int. rel.int."])
    for mz, intensity in peaks:
        relative = round(999 * intensity / maximum) if maximum else 0
        lines.append(f"  {mz:.15g} {intensity:.15g} {relative}")
    lines.append("//")
    return accession, "\n".join(lines) + "\n"


def convert_JSON_to_TXT(path_to_json_file, output_folder=None):
    """Export one TXT file per spectrum into a new folder; return its path."""
    output_folder = Path(output_folder) if output_folder is not None else (
        Path(__file__).parent / "txt_exports" / datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    )
    # Validate the input before creating the output folder.
    source_records = list(load_json_records(path_to_json_file))
    if not source_records:
        raise ValueError("Input contains no spectra")
    output_folder.mkdir(parents=True, exist_ok=False)

    for number, (_, record) in enumerate(source_records, start=1):
        accession, content = convert_record_to_txt(record, number)
        output_file = output_folder / f"{accession}.txt"
        output_file.write_text(content, encoding="utf-8")

    print(f"Exported {len(source_records)} TXT files to {output_folder}")
    return output_folder


def convert_record_to_sdf(record, number):
    """Build one RDKit molecule with MoNA-style spectrum properties.

    Use SMILES, falling back to InChI, and generate 2D coordinates.
    Preserve source intensities (no normalization) and all source metadata.
    A molecule is required: invalid/missing structures raise ValueError.
    Each source spectrum remains a separate record, even for the same compound.
    """
    from rdkit import Chem
    from rdkit.Chem import rdDepictor

    peaks = validated_peaks(record)
    molecule = None
    structure_source = None
    for field, parser in (("smiles", Chem.MolFromSmiles), ("inchi", Chem.MolFromInchi)):
        value = record.get(field)
        if isinstance(value, str) and value.strip():
            try:
                molecule = parser(value.strip())
            except (ValueError, RuntimeError):
                molecule = None
            if molecule is not None and molecule.GetNumAtoms() > 0:
                structure_source = field.upper()
                break
            molecule = None
    if molecule is None:
        raise ValueError("No valid molecular structure: supply a valid SMILES or InChI")

    rdDepictor.Compute2DCoords(molecule)
    name = text(record.get("compound_name") or "Unknown")
    molecule.SetProp("_Name", name)

    def add(label, value):
        if value is not None and value != "":
            rendered = text(value)
            # A standalone SDF delimiter cannot be used as a property value.
            if rendered.strip() == "$$$$":
                rendered = json.dumps(rendered)
            molecule.SetProp(label, rendered)

    add("NAME", name)
    fields = {
        "PRECURSOR TYPE": "adduct", "PRECURSOR M/Z": "precursor_mz",
        "INSTRUMENT TYPE": "instrument_type", "INSTRUMENT": "instrument",
        "COLLISION ENERGY": "collision_energy", "INCHIKEY": "inchikey",
        "FORMULA": "formula", "EXACT MASS": "exact_mass",
        "SMILES": "smiles", "INCHI": "inchi", "RETENTION TIME": "rt",
        "IONIZATION": "ion_source", "CHARGE": "charge", "MERGE TYPE": "merge_type",
        "CONTRIBUTOR": "investigator",
    }
    for label, field in fields.items():
        add(label, record.get(field))
    if record.get("ms_level") is not None:
        add("SPECTRUM TYPE", f"MS{record['ms_level']}")
    add("ION MODE", {"+": "P", "-": "N"}.get(record.get("polarity"), record.get("polarity")))
    add("ID", f"LOCAL{number:06d}")
    add("STRUCTURE SOURCE", structure_source)
    add("FRAGMENTATION METHOD", record.get("fragmentation_method") or record.get("fragmenation_method"))
    add("NUM PEAKS", len(peaks))
    # repr(float) preserves numeric precision when the SDF is read back.
    molecule.SetProp("MASS SPECTRAL PEAKS", "\n".join(
        f"{mz!r} {intensity!r}" for mz, intensity in peaks
    ))
    add("COMMENT", "Local export using MoNA-style field names. Intensities preserved from source. "
        "Collision-energy and retention-time units are unspecified in the source. "
        "SOURCE_METADATA retains the original metadata, including list-valued fields.")
    add("SOURCE_METADATA", {key: value for key, value in record.items() if key != "peaks"})
    return molecule


def convert_JSON_to_SDF(path_to_json_file, output_file=None):
    """Write all spectra into one SDF, without overwriting an existing file.

    Stream records to a temporary file, publishing the SDF only after success.
    Errors identify the input line and compound; no spectra are silently skipped.
    Return the completed SDF path.
    """
    try:
        from rdkit import Chem
    except ImportError as error:
        raise RuntimeError("SDF export requires RDKit: install it with python -m pip install rdkit") from error

    input_path = Path(path_to_json_file)
    output_file = Path(output_file) if output_file is not None else (
        Path(__file__).parent / "sdf_exports" /
        f"{input_path.stem}_{datetime.now():%Y%m%d_%H%M%S_%f}.sdf"
    )
    if output_file.exists():
        raise FileExistsError(f"Output already exists: {output_file}")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = None
    count = 0
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="\n",
                                         suffix=".tmp", dir=output_file.parent, delete=False) as temporary:
            temporary_path = Path(temporary.name)
            with Chem.SDWriter(temporary.file) as writer:
                for line_number, record in load_json_records(input_path):
                    try:
                        molecule = convert_record_to_sdf(record, count + 1)
                        writer.write(molecule)
                    except (ValueError, KeyError, TypeError, RuntimeError, OverflowError) as error:
                        name = text(record.get("compound_name") or "Unknown")
                        raise ValueError(f"Input line {line_number} ({name}): {error}") from error
                    count += 1
        if not count:
            raise ValueError("Input contains no spectra")
        # Linking the completed file fails if the destination already exists.
        # Keeping both paths on the same filesystem makes publication atomic.
        output_file.hardlink_to(temporary_path)
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
    print(f"Exported {count} spectra to {output_file}")
    return output_file


if __name__ == "__main__":
    # PyCharm Run uses this format and the sample input unless arguments are supplied.
    OUTPUT_FORMAT = "sdf"  # Change to "txt" for MassBank-style text export.
    parser = argparse.ArgumentParser(description="Convert MZmine JSON Lines to SDF or TXT.")
    parser.add_argument("input", nargs="?", type=Path,
                        default=Path(__file__).parent / "testing_files" / "testting_file.json")
    parser.add_argument("--format", choices=("sdf", "txt"), default=OUTPUT_FORMAT)
    parser.add_argument("--output", type=Path, help="SDF file path, or a new directory for TXT")
    args = parser.parse_args()
    try:
        if args.format == "sdf":
            convert_JSON_to_SDF(args.input, args.output)
        else:
            convert_JSON_to_TXT(args.input, args.output)
    except (ValueError, KeyError, TypeError, RuntimeError, OSError) as error:
        parser.exit(1, f"Conversion failed: {error}\n")
