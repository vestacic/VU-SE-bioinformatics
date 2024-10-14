from pathlib import Path

import numpy as np
from Bio import SeqIO

from organism import Organism


def get_input_data_files_paths(
    root_directory: Path, file_index_range: int
) -> list[Path]:
    """
    Constructs list of Path objects to load sequences from.
    """

    bacterial_files_paths = []
    mamalian_file_paths = []
    for i in range(1, file_index_range + 1):
        bacterial_files_paths.append(
            root_directory / f"input_data/bacterial{i}.fasta"
        )
        mamalian_file_paths.append(
            root_directory / f"input_data/mamalian{i}.fasta"
        )
    return bacterial_files_paths + mamalian_file_paths


def load_organism_from_fasta_file(file_path: Path) -> Organism:
    """
    Loads given .fasta file into Organism object.
    Expects 1 DNA sequence in file.
    """

    record = next(SeqIO.parse(file_path, "fasta"))
    return Organism(
        name=record.id, dna_sequence=np.array(record.seq, dtype="<U1")
    )


def export_distance_matrix_in_phylip(
    output_dir: Path,
    file_name: Path,
    matrix_size: int,
    distance_matrix: list[int, int],
    names: list[str],
) -> None:
    """
    Exports given distance matrix to .phy file.
    Creates output directory if it does not exist already.
    """

    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / (file_name + ".phy"), "w") as file:
        file.write(f"{matrix_size} ")

        for i in range(0, matrix_size):
            file.write(
                "\n" + names[i] + " " + " ".join(map(str, distance_matrix[i]))
            )
