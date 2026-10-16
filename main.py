
import json
from pathlib import Path
from datetime import date, datetime
import math

import numpy as np
from matchms import Spectrum

def text(value):
    if isinstance(value, (list, dict)):
        value = json.dumps(value, ensure_ascii=False)
    return str(value).replace("\r", " ").replace("\n", " ")


def convert_record(record, number):
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
    peaks = sorted((float(mz), float(intensity)) for mz, intensity in record["peaks"])
    if not peaks:
        raise ValueError("Spectrum has no peaks")
    if any(not math.isfinite(mz) or not math.isfinite(i) or mz <= 0 or i < 0
           for mz, i in peaks):
        raise ValueError("Peaks must have positive finite m/z and nonnegative finite intensity")
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


def convert_JSON_to_TXT(path_to_json_file):

    spectra = []
    source_records = []

    with open(path_to_json_file, encoding="utf-8-sig") as file:
        for line in file:
            if not line.strip():
                continue

            record = json.loads(line)
            source_record = record.copy()
            peaks = np.asarray(record.pop("peaks"), dtype=float)

            if peaks.size == 0:
                continue

            source_records.append(source_record)

            # Matchms requires peaks sorted by increasing m/z.
            peaks = peaks[np.argsort(peaks[:, 0])]

            spectra.append(
                Spectrum(
                    mz=peaks[:, 0],
                    intensities=peaks[:, 1],
                    metadata=record,
                )
            )

    # Export the loaded spectra to one MassBank-style TXT file per spectrum.
    # Each run gets a new folder, so earlier exports are kept.
    output_folder = Path(__file__).parent / "txt_exports" / datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    output_folder.mkdir(parents=True, exist_ok=False)

    for number, record in enumerate(source_records, start=1):
        accession, content = convert_record(record, number)
        output_file = output_folder / f"{accession}.txt"
        output_file.write_text(content, encoding="utf-8")

    print(f"Exported {len(source_records)} TXT files to {output_folder}")


if __name__ == "__main__":
    # Read the input file next to this script.
    file_path = Path(__file__).parent / "testing_files" / "testting_file.json"
    convert_JSON_to_TXT(file_path)
