import numpy as np
from Bio.Seq import Seq


class SequenceHandler:
    NUCLEOTIDES = np.array(["A", "C", "G", "T"], dtype="<U1")
    COMPLEMENTS = np.array(["T", "G", "C", "A"], dtype="<U1")
    START_CODON = np.array(["A", "T", "G"])
    STOP_CODONS = [
        np.array(["T", "A", "A"]),
        np.array(["T", "A", "G"]),
        np.array(["T", "G", "A"]),
    ]

    @classmethod
    def calculate_all_possible_sequence_frequencies(
        cls,
        dna_sequence: np.ndarray,
        codon_frequency: dict,
        dicodon_frequency: dict,
    ) -> None:
        """
        Calculates codon and dicodon frequencies
        for both sequence and its reverse complement
        in all reading frames.
        """

        cls.calculate_sequence_frequencies(
            dna_sequence=dna_sequence,
            codon_frequency=codon_frequency,
            dicodon_frequency=dicodon_frequency,
        )
        cls.calculate_sequence_frequencies(
            dna_sequence=dna_sequence[1:],
            codon_frequency=codon_frequency,
            dicodon_frequency=dicodon_frequency,
        )
        cls.calculate_sequence_frequencies(
            dna_sequence=dna_sequence[2:],
            codon_frequency=codon_frequency,
            dicodon_frequency=dicodon_frequency,
        )
        reverse_complement_dna_sequence = cls.calculate_reverse_complement(
            dna_sequence=dna_sequence
        )
        cls.calculate_sequence_frequencies(
            dna_sequence=reverse_complement_dna_sequence,
            codon_frequency=codon_frequency,
            dicodon_frequency=dicodon_frequency,
        )
        cls.calculate_sequence_frequencies(
            dna_sequence=reverse_complement_dna_sequence[1:],
            codon_frequency=codon_frequency,
            dicodon_frequency=dicodon_frequency,
        )
        cls.calculate_sequence_frequencies(
            dna_sequence=reverse_complement_dna_sequence[2:],
            codon_frequency=codon_frequency,
            dicodon_frequency=dicodon_frequency,
        )

    @classmethod
    def calculate_sequence_frequencies(
        cls,
        dna_sequence: np.ndarray,
        codon_frequency: dict,
        dicodon_frequency: dict,
    ) -> None:
        """
        Calculates codon and dicodon frequencies for a given sequence.
        """

        all_start_stop_codon_pairs = cls.find_start_stop_codon_pairs(
            dna_sequence=dna_sequence
        )
        longest_start_stop_codon_pairs_for_each_stop_codon = (
            cls.filter_longest_pair_for_each_codon(
                start_stop_codon_pairs=all_start_stop_codon_pairs
            )
        )
        start_stop_codon_pairs_not_smaller_than_100 = cls.filter_pairs_smaller_than(
            start_stop_codon_pairs=longest_start_stop_codon_pairs_for_each_stop_codon
        )
        cut_dna_sequences = cls.cut_sequences(
            initial_dna_sequence=dna_sequence,
            start_stop_codon_pairs=start_stop_codon_pairs_not_smaller_than_100,
        )
        protein_sequences = cls.translate_sequences(
            dna_sequences=cut_dna_sequences
        )
        cls.calculate_codon_frequence(
            codon_frequency=codon_frequency,
            protein_sequences=protein_sequences,
        )
        cls.calculate_dicodon_frequence(
            dicodon_frequency=dicodon_frequency,
            protein_sequences=protein_sequences,
        )

    @classmethod
    def calculate_reverse_complement(
        cls, dna_sequence: np.ndarray
    ) -> np.ndarray:
        """
        Finds reverse complement for given sequence.
        """

        index_map = np.searchsorted(cls.NUCLEOTIDES, dna_sequence)
        reverse_complement = cls.COMPLEMENTS[index_map]
        return reverse_complement[::-1]

    @classmethod
    def find_start_stop_codon_pairs(
        cls, dna_sequence: np.ndarray
    ) -> list[tuple[int, int]]:
        """
        Finds all start and stop codon pairs that do not have stop codon in between.
        """

        sequence_length = len(dna_sequence)

        pairs = []
        found_start_codons = []

        for i in range(0, sequence_length, 3):
            codon = dna_sequence[i : i + 3]

            if np.array_equal(codon, cls.START_CODON):
                found_start_codons.append(i)

            elif any(
                np.array_equal(codon, stop_codon)
                for stop_codon in cls.STOP_CODONS
            ):
                pairs.extend(
                    (start_index, i) for start_index in found_start_codons
                )
                found_start_codons.clear()

        return np.array(pairs, dtype=int)

    @classmethod
    def filter_longest_pair_for_each_codon(
        cls, start_stop_codon_pairs: np.ndarray
    ) -> dict:
        """
        For each stop codon filters out the start codon
        that is furthest from it and does not have another stop codon in between.
        """

        starts = start_stop_codon_pairs[:, 0]
        stops = start_stop_codon_pairs[:, 1]

        unique_stops, indices = np.unique(stops, return_inverse=True)

        longest_starts = np.full(len(unique_stops), -1, dtype=int)

        for i, stop in enumerate(unique_stops):
            relevant_starts = starts[stops == stop]
            longest_starts[i] = np.min(relevant_starts)

        return np.column_stack((longest_starts, unique_stops))

    @classmethod
    def filter_pairs_smaller_than(
        cls, start_stop_codon_pairs: np.ndarray, min_length: int = 100
    ) -> dict:
        """
        Filters out start and stop codon pairs that have length less than 100 bp.
        """

        lengths = (
            start_stop_codon_pairs[:, 1] - start_stop_codon_pairs[:, 0] + 3
        )

        filtered_pairs = start_stop_codon_pairs[lengths >= min_length]

        return filtered_pairs

    @classmethod
    def cut_sequences(
        cls,
        initial_dna_sequence: np.ndarray,
        start_stop_codon_pairs: np.ndarray,
    ) -> np.ndarray:
        """
        Cuts sequences from given indices.
        """

        cut_dna_sequences = np.empty(len(start_stop_codon_pairs), dtype=object)
        for i, (start, stop) in enumerate(start_stop_codon_pairs):
            cut_dna_sequences[i] = "".join(
                initial_dna_sequence[start : stop + 3]
            )

        return cut_dna_sequences

    @classmethod
    def translate_sequences(cls, dna_sequences: np.ndarray) -> np.ndarray:
        protein_sequences = np.empty(len(dna_sequences), dtype=object)
        for i, sequence in enumerate(dna_sequences):
            dna_sequence = Seq(sequence)
            protein_sequence = str(dna_sequence.translate())
            protein_sequences[i] = protein_sequence

        return protein_sequences

    @classmethod
    def calculate_codon_frequence(
        cls, codon_frequency: dict, protein_sequences: np.ndarray
    ) -> None:
        for protein_sequence in protein_sequences:
            for codon in protein_sequence[:-1]:
                codon_frequency[codon] += 1

    @classmethod
    def calculate_dicodon_frequence(
        cls, dicodon_frequency: dict, protein_sequences: np.ndarray
    ) -> None:
        for protein_sequence in protein_sequences:
            for i in range(len(protein_sequence) - 2):
                dicodon = protein_sequence[i] + protein_sequence[i + 1]
                dicodon_frequency[dicodon] += 1
