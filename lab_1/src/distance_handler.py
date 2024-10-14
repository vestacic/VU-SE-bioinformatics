from decimal import ROUND_HALF_UP, Decimal


def calculate_chi_square_distance(
    frequencies_1: dict, frequencies_2: dict
) -> Decimal:
    """
    Calculates distance between 2 given codon frequencies.
    Rounds up to 4 digits after floating point.
    """

    chi_square_distance = Decimal(0)
    for codon in frequencies_1:
        frequency_1 = frequencies_1[codon]
        frequency_2 = frequencies_2[codon]
        if frequency_1 + frequency_2 != 0:
            chi_square_distance += ((frequency_1 - frequency_2) ** 2) / (
                frequency_1 + frequency_2
            )

    return chi_square_distance.quantize(
        Decimal("0.0001"), rounding=ROUND_HALF_UP
    )
