from pathlib import Path

import numpy as np
from Bio import SeqIO

from organism import Organism


def load_organism_from_fasta_file(file_path: Path) -> Organism:
    record = next(SeqIO.parse(file_path, "fasta"))
    return Organism(
        name=record.id, dna_sequence=np.array(record.seq, dtype="<U1")
    )
