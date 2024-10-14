import itertools
from decimal import Decimal

import numpy as np

from sequence_handler import SequenceHandler


class Organism:
    AMINO_ACIDS = [
        "A",
        "R",
        "N",
        "D",
        "C",
        "E",
        "Q",
        "G",
        "H",
        "I",
        "L",
        "K",
        "M",
        "F",
        "P",
        "S",
        "T",
        "W",
        "Y",
        "V",
    ]

    def __init__(self, name: str, dna_sequence: np.ndarray) -> None:
        self.name = name
        self.dna_sequence = dna_sequence

        self.amino_acid_codon_frequency = (
            Organism.create_empty_codon_frequency_dict()
        )

        self.amino_acid_dicodon_frequency = (
            Organism.create_empty_dicodon_frequency_dict()
        )

        self.amino_acid_absolute_codon_frequency = (
            Organism.create_empty_codon_frequency_dict()
        )

        self.amino_acid_absolute_dicodon_frequency = (
            Organism.create_empty_dicodon_frequency_dict()
        )

    @classmethod
    def create_empty_dicodon_frequency_dict(cls) -> dict:
        return {
            amino_acid_1 + amino_acid_2: 0
            for amino_acid_1, amino_acid_2 in itertools.product(
                cls.AMINO_ACIDS, cls.AMINO_ACIDS
            )
        }

    @classmethod
    def create_empty_codon_frequency_dict(cls) -> dict:
        return {amino_acid: 0 for amino_acid in cls.AMINO_ACIDS}

    def calculate_amino_acid_frequency(self) -> None:
        if self.amino_acid_codon_frequency:
            self.amino_acid_codon_frequency = (
                Organism.create_empty_codon_frequency_dict()
            )
        if self.amino_acid_dicodon_frequency:
            self.amino_acid_dicodon_frequency = (
                Organism.create_empty_dicodon_frequency_dict()
            )

        SequenceHandler.calculate_all_possible_sequence_frequencies(
            dna_sequence=self.dna_sequence,
            codon_frequency=self.amino_acid_codon_frequency,
            dicodon_frequency=self.amino_acid_dicodon_frequency,
        )

    def calculate_normalized_amino_acid_frequency(self) -> None:
        self.calculate_amino_acid_frequency()

        total_codon_number = sum(self.amino_acid_codon_frequency.values())
        self.amino_acid_absolute_codon_frequency = {
            codon: Decimal(frequency) / Decimal(total_codon_number)
            for codon, frequency in self.amino_acid_codon_frequency.items()
        }

        total_dicodon_number = sum(self.amino_acid_dicodon_frequency.values())
        self.amino_acid_absolute_dicodon_frequency = {
            codon: Decimal(frequency) / Decimal(total_dicodon_number)
            for codon, frequency in self.amino_acid_dicodon_frequency.items()
        }
