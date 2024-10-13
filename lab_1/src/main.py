from pathlib import Path

import file_handler
from distance_handler import calculate_chi_square_distance

DATA_FILES_INDEX_RANGE = 4


def main() -> None:
    root_directory_path = Path(__file__).parents[1]
    data_files_paths = file_handler.get_input_data_files_paths(
        root_directory=root_directory_path,
        file_index_range=DATA_FILES_INDEX_RANGE,
    )
    organisms = [
        file_handler.load_organism_from_fasta_file(file_path=data_file_path)
        for data_file_path in data_files_paths
    ]

    [
        organism.calculate_absolute_amino_acid_frequency()
        for organism in organisms
    ]

    codon_distance_matrix = [
        [0 for i in range(DATA_FILES_INDEX_RANGE * 2)]
        for j in range(DATA_FILES_INDEX_RANGE * 2)
    ]

    dicodon_distance_matrix = [
        [0 for i in range(DATA_FILES_INDEX_RANGE * 2)]
        for j in range(DATA_FILES_INDEX_RANGE * 2)
    ]

    for i in range(DATA_FILES_INDEX_RANGE * 2):
        for j in range(i, DATA_FILES_INDEX_RANGE * 2):
            codon_chi_square_distance = calculate_chi_square_distance(
                frequencies_1=organisms[i].amino_acid_absolute_codon_frequency,
                frequencies_2=organisms[j].amino_acid_absolute_codon_frequency,
            )

            dicodon_chi_square_distance = calculate_chi_square_distance(
                frequencies_1=organisms[
                    i
                ].amino_acid_absolute_dicodon_frequency,
                frequencies_2=organisms[
                    j
                ].amino_acid_absolute_dicodon_frequency,
            )

            codon_distance_matrix[i][j] = codon_distance_matrix[j][i] = (
                codon_chi_square_distance
            )
            dicodon_distance_matrix[i][j] = dicodon_distance_matrix[j][i] = (
                dicodon_chi_square_distance
            )
    organism_names = [organism.name for organism in organisms]

    output_dir = root_directory_path / "output_data"
    file_handler.export_distance_matrix_in_phylip(
        output_dir=output_dir,
        file_name="codon_distance_matrix.phylip",
        matrix_size=DATA_FILES_INDEX_RANGE * 2,
        distance_matrix=codon_distance_matrix,
        names=organism_names,
    )

    file_handler.export_distance_matrix_in_phylip(
        output_dir=output_dir,
        file_name="dicodon_distance_matrix.phylip",
        matrix_size=DATA_FILES_INDEX_RANGE * 2,
        distance_matrix=dicodon_distance_matrix,
        names=organism_names,
    )


if __name__ == "__main__":
    main()
