from decimal import Decimal


def calculate_chi_square_distance(
    frequencies_1: dict, frequencies_2: dict
) -> Decimal:
    chi_square_distance = Decimal(0)
    for codon in frequencies_1:
        frequency_1 = frequencies_1[codon]
        frequency_2 = frequencies_2[codon]
        if frequency_1 + frequency_2 != 0:
            chi_square_distance += ((frequency_1 - frequency_2) ** 2) / (
                frequency_1 + frequency_2
            )

    return chi_square_distance
