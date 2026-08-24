import argparse
from pathlib import Path
import re
import sys

from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import rdBase


SMI_HEADER_TOKENS = {"smiles", "smi"}


def is_smi_header(line: str) -> bool:
    columns = line.split()
    return bool(columns) and columns[0].casefold() in SMI_HEADER_TOKENS


def read_smi(path: Path) -> tuple[list[Chem.Mol], list[str]]:
    molecules = []
    errors = []
    first_record_checked = False
    record_index = 0

    with path.open(encoding="utf-8-sig") as file:
        for line_number, line in enumerate(file, start=1):
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            if not first_record_checked:
                first_record_checked = True

                if is_smi_header(line):
                    continue

            record_index += 1
            columns = line.split(maxsplit=1)
            smiles = columns[0]
            name = columns[1].strip() if len(columns) == 2 else ""

            if not name:
                name = f"ligand_{record_index:04d}"

            try:
                with rdBase.BlockLogs():
                    molecule = Chem.MolFromSmiles(smiles)
            except Exception as error:
                errors.append(
                    f"line {line_number} ({name}): "
                    f"{type(error).__name__}: {error}"
                )
                continue

            if molecule is None:
                errors.append(
                    f"line {line_number} ({name}): cannot parse SMILES {smiles!r}"
                )
                continue

            molecule.SetProp("_Name", name)
            molecules.append(molecule)

    return molecules, errors


def read_sdf(path: Path) -> tuple[list[Chem.Mol], list[str]]:
    supplier = Chem.SDMolSupplier(
        str(path),
        removeHs=False,
        sanitize=True,
    )

    molecules = []
    errors = []

    for index, molecule in enumerate(supplier, start=1):
        if molecule is None:
            errors.append(f"record {index}: cannot parse SDF molecule")
            continue

        name = ""

        if molecule.HasProp("_Name"):
            name = molecule.GetProp("_Name").strip()

        if not name and molecule.HasProp("name"):
            name = molecule.GetProp("name").strip()

        if not name:
            name = f"ligand_{index:04d}"

        molecule.SetProp("_Name", name)
        molecules.append(molecule)

    return molecules, errors


def read_molecules(path: Path) -> tuple[list[Chem.Mol], list[str]]:
    suffix = path.suffix.lower()

    if suffix == ".smi":
        return read_smi(path)

    if suffix == ".sdf":
        return read_sdf(path)

    raise ValueError(
        f"Unsupported input format: {suffix}. Supported formats: .smi, .sdf"
    )


def has_3d_coordinates(molecule: Chem.Mol) -> bool:
    if molecule.GetNumConformers() == 0:
        return False

    conformer = molecule.GetConformer()

    if conformer.Is3D():
        return True

    z_coordinates = [
        conformer.GetAtomPosition(index).z
        for index in range(molecule.GetNumAtoms())
    ]

    return bool(z_coordinates) and max(z_coordinates) - min(z_coordinates) > 1e-3


def generate_3d(molecule: Chem.Mol) -> Chem.Mol:
    molecule_copy = Chem.Mol(molecule)

    if has_3d_coordinates(molecule_copy):
        return Chem.AddHs(molecule_copy, addCoords=True)

    molecule_3d = Chem.AddHs(molecule_copy)
    molecule_3d.RemoveAllConformers()

    parameters = AllChem.ETKDGv3()

    with rdBase.BlockLogs():
        embedding_status = AllChem.EmbedMolecule(
            molecule_3d,
            parameters,
        )

    if embedding_status != 0:
        raise RuntimeError(f"3D embedding failed for {molecule.GetProp('_Name')}")

    return molecule_3d


def convert_to_pdbqt(molecule: Chem.Mol) -> str:
    preparator = MoleculePreparation(
        rigid_macrocycles=False,
        min_ring_size=6,
    )
    molecule_setups = preparator.prepare(molecule)

    if not molecule_setups:
        raise RuntimeError(f"Meeko preparation failed for {molecule.GetProp('_Name')}")

    pdbqt_string, success, error_message = PDBQTWriterLegacy.write_string(
        molecule_setups[0]
    )

    if not success:
        raise RuntimeError(
            f"PDBQT writing failed for {molecule.GetProp('_Name')}: {error_message}"
        )

    return pdbqt_string


def sanitize_filename(name: str) -> str:
    name = re.sub(r'[<>:"/\\|?*\x00-\x1f]', "_", name)
    name = name.strip().rstrip(". ")
    return name or "unnamed_ligand"


def make_output_paths(
    molecules: list[Chem.Mol],
    output_dir: Path,
) -> list[Path]:
    output_paths = []
    used_filenames = set()

    for molecule in molecules:
        base_name = sanitize_filename(molecule.GetProp("_Name"))
        candidate = base_name
        suffix = 2

        while f"{candidate}.pdbqt".casefold() in used_filenames:
            candidate = f"{base_name}_{suffix}"
            suffix += 1

        filename = f"{candidate}.pdbqt"
        used_filenames.add(filename.casefold())
        output_paths.append(output_dir / filename)

    return output_paths


def process_molecules(
    molecules: list[Chem.Mol],
    output_dir: Path,
) -> tuple[int, int]:
    output_dir.mkdir(parents=True, exist_ok=True)
    output_paths = make_output_paths(molecules, output_dir)

    successful = 0
    failed = 0
    total = len(molecules)

    for index, (molecule, output_path) in enumerate(
        zip(molecules, output_paths),
        start=1,
    ):
        name = molecule.GetProp("_Name")

        try:
            molecule_3d = generate_3d(molecule)
            pdbqt_string = convert_to_pdbqt(molecule_3d)

            output_path.write_text(pdbqt_string, encoding="utf-8")

            successful += 1
            print(f"[{index}/{total}] OK: {name} -> {output_path.name}")

        except Exception as error:
            failed += 1
            print(
                f"[{index}/{total}] FAILED: {name}: "
                f"{type(error).__name__}: {error}",
                file=sys.stderr,
            )

    return successful, failed


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare multiple ligands from SMI or SDF "
            "for AutoDock Vina in PDBQT format."
        )
    )

    parser.add_argument(
        "input",
        type=Path,
        help="Path to an input .smi or .sdf file.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for generated PDBQT files.",
    )

    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()

    input_path = arguments.input.resolve()
    output_dir = arguments.output_dir.resolve()

    if not input_path.is_file():
        print(f"ERROR: input file does not exist: {input_path}", file=sys.stderr)
        return

    try:
        molecules, read_errors = read_molecules(input_path)
    except (OSError, ValueError) as error:
        print(f"ERROR: {type(error).__name__}: {error}", file=sys.stderr)
        return

    if not molecules and not read_errors:
        read_errors.append("no molecule records found")

    for error in read_errors:
        print(f"READ FAILED: {error}", file=sys.stderr)

    successful, processing_failed = process_molecules(
        molecules,
        output_dir,
    )

    failed = len(read_errors) + processing_failed

    print()
    print(f"Input file: {input_path}")
    print(f"Output directory: {output_dir}")
    print(f"Successfully processed: {successful}")
    print(f"Failed: {failed}")


if __name__ == "__main__":
    main()
