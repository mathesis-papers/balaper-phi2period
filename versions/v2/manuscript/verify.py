import argparse
from functools import lru_cache
from fractions import Fraction
import json
from math import comb, factorial, lcm
from pathlib import Path


def trim(coefficients):
    while len(coefficients) > 1 and coefficients[-1] == 0:
        coefficients.pop()
    return coefficients


def add(left, right):
    return trim([
        (left[index] if index < len(left) else 0)
        + (right[index] if index < len(right) else 0)
        for index in range(max(len(left), len(right)))
    ])


def scale(coefficients, scalar):
    return trim([scalar * value for value in coefficients])


def multiply(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            result[left_index + right_index] += left_value * right_value
    return trim(result)


def power(coefficients, exponent):
    result = [1]
    for _ in range(exponent):
        result = multiply(result, coefficients)
    return result


def derivative(coefficients):
    return trim([index * coefficients[index] for index in range(1, len(coefficients))] or [0])


def column_operator(coefficients, parameter, parity):
    intermediate = add(scale(coefficients, parity), scale(multiply([1, 1], derivative(coefficients)), 2))
    first = multiply([1, parity + 2 * parameter - 1], intermediate)
    second = scale(multiply([0, 1, 1], derivative(intermediate)), 2)
    return multiply(power([1, 1], 2 * parameter - 2), add(first, second))


def identity(dimension):
    return [[int(row == column) for column in range(dimension)] for row in range(dimension)]


def matrix_multiply(left, right, modulus):
    dimension = len(left)
    result = [[0] * dimension for _ in range(dimension)]
    for row in range(dimension):
        for intermediate, value in enumerate(left[row]):
            if value:
                for column in range(dimension):
                    result[row][column] += value * right[intermediate][column]
        result[row] = [value % modulus for value in result[row]]
    return result


def matrix_power(matrix, exponent, modulus):
    result = identity(len(matrix))
    while exponent:
        if exponent & 1:
            result = matrix_multiply(result, matrix, modulus)
        matrix = matrix_multiply(matrix, matrix, modulus)
        exponent //= 2
    return result


def apply(matrix, vector, modulus):
    return [sum(value * coefficient for value, coefficient in zip(row, vector)) % modulus for row in matrix]


def column_matrix(parameter, parity, modulus):
    dimension = (modulus + 1) // 2
    columns = []
    for column in range(dimension):
        coefficients = column_operator([0] * column + [1], parameter, parity)
        columns.append([(coefficients[row] if row < len(coefficients) else 0) % modulus for row in range(dimension)])
    return [[columns[column][row] for column in range(dimension)] for row in range(dimension)]


def column_residues(parameter, exponent, modulus, count):
    matrix = column_matrix(parameter, exponent % 2, modulus)
    coefficients = power([1, 1], exponent // 2)
    state = [(coefficients[index] if index < len(coefficients) else 0) % modulus for index in range(len(matrix))]
    values = []
    for _ in range(count):
        values.append(state[0])
        state = apply(matrix, state, modulus)
    return values


@lru_cache(maxsize=None)
def recurrence(parameter, exponent, index):
    if index == 0:
        return 1
    return exponent * (
        (exponent + 2 * parameter - 1) * recurrence(parameter, exponent + 4 * parameter - 2, index - 1)
        - (exponent + 2 * parameter - 2) * recurrence(parameter, exponent + 4 * parameter - 4, index - 1)
    )


def coefficient_checks():
    families = [
        ('A281180', 2, 4, 1, [1, 4, 88, 4672, 454144, 70084096, 15728822272, 4836914249728]),
        ('A281181', 2, 1, 0, [1, 1, 13, 493, 37369, 4732249, 901188997, 240798388357]),
        ('A281183', 2, 2, 0, [1, 2, 32, 1376, 114176, 15519488, 3132551168, 879422726144]),
        ('A281184', 2, 3, 0, [1, 3, 57, 2739, 246801, 35822307, 7636142793, 2246286827091]),
        ('A281440', 4, 1, 0, [1, 1, 29, 2605, 478745, 148838425, 70130095925, 46612385237125]),
    ]
    results = []
    for sequence, parameter, exponent, offset, expected in families:
        polynomial = power([1, 1], exponent // 2)
        values = []
        for index in range(21):
            value = polynomial[0]
            assert value == recurrence(parameter, exponent, index)
            values.append(value)
            polynomial = column_operator(polynomial, parameter, exponent % 2)
        assert values[:len(expected)] == expected
        results.append({'sequence': sequence, 'first_index': offset, 'last_index': offset + 20,
                        'methods': ['quadratic derivation in s^2', 'fixed-column integer recurrence'],
                        'oeis_prefix_length': len(expected), 'prefix': expected})
    coefficients = [0]
    for degree in range(40):
        square = [sum(comb(index, part) * coefficients[part] * coefficients[index - part]
                      for part in range(index + 1)) for index in range(degree + 1)]
        fourth = sum(comb(degree, part) * square[part] * square[degree - part] for part in range(degree + 1))
        coefficients.append(int(degree == 0) + 2 * square[degree] + fourth)
    polynomial = [0, 1]
    for degree in range(41):
        assert polynomial[0] == coefficients[degree]
        polynomial = multiply([1, 0, 2, 0, 1], derivative(polynomial))
    for index in range(1, 21):
        assert coefficients[2 * index - 1] == recurrence(2, 4, index - 1)
    return {'families': results, 'main_ode_and_polynomial_derivation': {'derivative_indices': [0, 40], 'equal': True}}


def prime_divisors(integer):
    return [candidate for candidate in range(2, integer + 1)
            if integer % candidate == 0 and all(candidate % divisor for divisor in range(2, candidate))]


def exact_period_certificates():
    specifications = [
        ('A281180', 2, 4, 7, 18, 1, 1),
        ('A281181', 2, 1, 10, 8, 0, 0),
        ('A281181', 2, 1, 9, 3, 1, 0),
        ('A281181', 2, 1, 13, 72, 0, 0),
        ('A281183', 2, 2, 11, 25, 1, 0),
        ('A281184', 2, 3, 11, 25, 1, 0),
        ('A281440', 4, 1, 13, 72, 0, 0),
    ]
    results = []
    for sequence, parameter, exponent, modulus, period, origin, offset in specifications:
        dimension = (modulus + 1) // 2
        start = origin - offset
        values = column_residues(parameter, exponent, modulus, start + period + dimension)
        differences = [(values[index + period] - values[index]) % modulus for index in range(start, start + dimension)]
        assert differences == [0] * dimension
        witnesses = []
        for prime in prime_divisors(period):
            smaller = period // prime
            witness = next(index for index in range(start, start + dimension)
                           if values[index + smaller] != values[index])
            witnesses.append({'excluded_shift': smaller, 'index': witness + offset,
                              'value': values[witness], 'shifted_value': values[witness + smaller]})
        results.append({'sequence': sequence, 'modulus': modulus, 'least_period': period,
                        'first_index': origin, 'dimension': dimension, 'matrix': column_matrix(parameter, exponent % 2, modulus),
                        'zero_observation_differences': differences, 'minimality_witnesses': witnesses,
                        'all_index_implication': 'Cayley-Hamilton for the displayed finite matrix'})
    values = column_residues(2, 1, 12, 7)
    assert values == [1] * 7
    return {'periods': results, 'A281181_mod_12': {'dimension': 6, 'values_at_0_through_6': values,
            'all_index_implication': 'Cayley-Hamilton applied to six successive zero differences'}}


def polynomial_operator_checks():
    polynomials = [[1], [1, 1], [1, 0, 2, 0, 1], [1, -2, 0, 3], [1, 0, 0, 0, 1]]
    prime_powers = [(2, 1), (3, 1), (2, 2), (5, 1), (7, 1), (2, 3),
                    (3, 2), (11, 1), (13, 1), (2, 4), (5, 2), (3, 3)]
    results = []
    for prime, exponent in prime_powers:
        modulus = prime ** exponent
        tail = prime ** (2 * exponent - 1)
        shift = ((prime - 1) * prime ** (exponent - 1)) ** 2
        for polynomial in polynomials:
            matrix = [[0] * modulus for _ in range(modulus)]
            for column in range(1, modulus):
                for degree, coefficient in enumerate(polynomial):
                    row = column - 1 + degree
                    if row < modulus:
                        matrix[row][column] = column * coefficient % modulus
            delayed = matrix_power(matrix, tail, modulus)
            shift_matrix = matrix_power(matrix, shift, modulus)
            difference = [[(shift_matrix[row][column] - int(row == column)) % modulus
                           for column in range(modulus)] for row in range(modulus)]
            product = matrix_multiply(delayed, difference, modulus)
            assert all(value == 0 for row in product for value in row)
            results.append({'prime': prime, 'exponent': exponent, 'modulus': modulus,
                            'polynomial_coefficients': polynomial, 'tail_derivative_index': tail,
                            'derivative_shift': shift, 'zero_matrix': True})
    return {'cases': len(results), 'identities': results, 'scope': 'Finite matrix checks of the polynomial-flow identity'}


def derivation_matrix(exponent, modulus):
    dimension = 2 * modulus
    matrix = [[0] * dimension for _ in range(dimension)]
    for parity in range(2):
        for degree in range(modulus):
            terms = []
            if degree:
                terms.append((degree, degree - 1, exponent + parity))
            if parity:
                terms.append((1, degree + 1, exponent - 1))
            for scalar, power_index, c_power in terms:
                for index in range(c_power // 2 + 1):
                    row_degree = power_index + 2 * index
                    if row_degree < modulus:
                        row = (c_power % 2) * modulus + row_degree
                        matrix[row][parity * modulus + degree] += scalar * comb(c_power // 2, index)
    return [[value % modulus for value in row] for row in matrix]


def subtract_scalar(matrix, scalar, modulus):
    return [[(value - scalar * int(row == column)) % modulus
             for column, value in enumerate(values)] for row, values in enumerate(matrix)]


def quadratic_operator_checks():
    results = []
    for exponent in range(1, 10):
        for prime, height in [(2, 1), (2, 2), (2, 3), (2, 4), (3, 1), (3, 2),
                              (3, 3), (5, 1), (5, 2), (7, 1), (11, 1), (13, 1)]:
            modulus = prime ** height
            matrix = derivation_matrix(exponent, modulus)
            square = matrix_multiply(matrix, matrix, modulus)
            totient = (prime - 1) * prime ** (height - 1)
            period = totient ** 2 // 2 if modulus >= 3 else 1
            onset = 2 ** (height - 1) if prime == 2 else (prime ** (2 * height - 1) + 1) // 2
            difference = subtract_scalar(matrix_power(square, period, modulus), 1, modulus)
            product = matrix_multiply(matrix_power(square, onset, modulus), difference, modulus)
            assert all(value == 0 for row in product for value in row)
            if prime == 2:
                row_width = 2 ** (height - 2) if height >= 2 else 1
                multiplier = 1
            else:
                row_width = (prime - 1) * prime ** (height - 1) // 2
                multiplier = (-1) ** ((prime - 1) // 2) % prime
                for index in range((prime - 1) // 2):
                    multiplier = multiplier * (exponent + 2 * index) * pow(2 * (index + 1), -1, prime) % prime
            row_difference = subtract_scalar(matrix_power(square, row_width, modulus), multiplier, modulus)
            row_product = matrix_multiply(matrix_power(row_difference, height, modulus),
                                          matrix_power(square, height, modulus), modulus)
            assert all(value == 0 for row in row_product for value in row)
            results.append({'exponent': exponent, 'modulus': modulus, 'rank': 2 * modulus,
                            'column_period': period, 'column_onset': onset, 'period_identity_zero': True,
                            'row_width': row_width, 'row_multiplier_lift': multiplier,
                            'row_order_and_onset': height, 'row_annihilator_zero': True})
    return {'cases': len(results), 'identities': results,
            'scope': 'Exact operator identities on the full finite quadratic algebra'}


def general_column_operator(coefficients, exponent, parity):
    intermediate = add(scale(coefficients, parity), scale(multiply([1, 1], derivative(coefficients)), 2))
    first = multiply([1, parity + exponent - 1], intermediate)
    second = scale(multiply([0, 1, 1], derivative(intermediate)), 2)
    return multiply(power([1, 1], exponent - 2), add(first, second))


@lru_cache(maxsize=None)
def general_recurrence(exponent, column, index):
    if index == 0:
        return 1
    if column == 0:
        return 0
    return column * (
        (column + exponent - 1) * general_recurrence(exponent, column + 2 * exponent - 2, index - 1)
        - (column + exponent - 2) * general_recurrence(exponent, column + 2 * exponent - 4, index - 1)
    )


def neighboring_coefficients():
    expected = [
        ('A281427', 5, 5, 1, [1, 5, 145, 10325, 1357825, 284963525, 87274812625, 36716097543125]),
        ('A281428', 5, 1, 0, [1, 1, 17, 865, 88865, 15335425, 3993275825, 1462392957025]),
        ('A281432', 6, 1, 0, [1, 1, 21, 1341, 173961, 38032281, 12572222301, 5853055718421]),
        ('A281435', 7, 7, 1, [1, 7, 301, 32347, 6476281, 2080072687, 978357441061, 633671918506627]),
        ('A281436', 7, 1, 0, [1, 1, 25, 1921, 301105, 79715041, 31953352585, 18055521444961]),
        ('A281439', 8, 8, 1, [1, 8, 400, 50240, 11785600, 4441856000, 2454008012800, 1868290095104000]),
        ('A281443', 9, 9, 1, [1, 9, 513, 73737, 19832769, 8579840841, 5445003346497, 4764370917991113]),
        ('A281444', 9, 1, 0, [1, 1, 33, 3393, 715329, 255419649, 138324673377, 105723078203457]),
    ]
    results = []
    for sequence, exponent, column, offset, prefix in expected:
        polynomial = power([1, 1], column // 2)
        values = []
        for index in range(13):
            assert polynomial[0] == general_recurrence(exponent, column, index)
            values.append(polynomial[0])
            polynomial = general_column_operator(polynomial, exponent, column % 2)
        assert values[:8] == prefix
        results.append({'sequence': sequence, 'exponent': exponent, 'column': column,
                        'first_index': offset, 'last_index': offset + 12,
                        'oeis_prefix_length': 8, 'prefix': prefix, 'equal': True})
    return results


def residue_row_certificates():
    specifications = [('A281180', 2, 4, 7, 6, 2, 1, 1),
                      ('A281181', 2, 1, 13, 6, 7, 0, 0),
                      ('A281183', 2, 2, 11, 5, 5, 1, 0),
                      ('A281184', 2, 3, 11, 5, 5, 1, 0),
                      ('A281440', 4, 1, 13, 6, 6, 1, 0)]
    results = []
    for sequence, parameter, column, modulus, width, multiplier, origin, offset in specifications:
        dimension = (modulus + 1) // 2
        start = origin - offset
        values = column_residues(parameter, column, modulus, start + width + dimension)
        differences = [(values[index + width] - multiplier * values[index]) % modulus
                       for index in range(start, start + dimension)]
        assert differences == [0] * dimension
        results.append({'sequence': sequence, 'modulus': modulus, 'row_width': width,
                        'row_multiplier': multiplier, 'first_index': origin,
                        'dimension': dimension, 'zero_observation_differences': differences,
                        'all_index_implication': 'Cayley-Hamilton for the column matrix'})
    tail = column_residues(2, 3, 9, 7)
    assert tail[1:] == [3] * 6
    return {'geometric_rows': results, 'A281184_mod_9': {'dimension': 5, 'first_index': 1,
            'values_at_1_through_6': tail[1:], 'all_index_implication': 'Cayley-Hamilton for five zero successive differences'}}


def divide_monic(dividend, divisor):
    dividend = dividend[:]
    quotient = [0] * (len(dividend) - len(divisor) + 1)
    for degree in range(len(quotient) - 1, -1, -1):
        quotient[degree] = dividend[degree + len(divisor) - 1]
        for index, value in enumerate(divisor):
            dividend[degree + index] -= quotient[degree] * value
    assert all(value == 0 for value in dividend)
    return trim(quotient)


def triangle_checks():
    polynomial = [1, 0, 2, 0, 1]
    previous = []
    rows = []
    for index in range(1, 13):
        quotient = divide_monic(polynomial, power([1, 0, 1], 2 * index))
        assert all(value == 0 for value in quotient[1::2])
        coefficients = quotient[::2]
        row = [1] if index == 1 else [
            2 * (column + 1) * (2 * column + 1) * (previous[column + 1] if column + 1 < len(previous) else 0)
            + 4 * ((4 * column + 1) * index + 2 * column ** 2 - 3 * column - 1) * (previous[column] if column < len(previous) else 0)
            + 2 * (4 * index + 2 * column - 3) * (2 * index + column - 3) * (previous[column - 1] if column else 0)
            for column in range(index)]
        assert coefficients == row
        current = recurrence(2, 4, index - 1)
        following = recurrence(2, 4, index)
        second = (following - 4 * index * current) // 2
        assert following - 4 * index * current == 2 * second
        assert row[0] == current and (row[1] if len(row) > 1 else 0) == second
        rows.append({'index': index, 'row': row, 'second_column': second})
        previous = row
        for _ in range(2):
            polynomial = multiply([1, 0, 2, 0, 1], derivative(polynomial))
    return {'indices': [1, 12], 'derivative_polynomial_equals_row_recurrence': True, 'rows': rows}


def derivative_formula_checks():
    maximum = 12
    cosine = [Fraction((-1) ** (degree // 2), factorial(degree)) if degree % 2 == 0 else Fraction(0)
              for degree in range(maximum + 1)]
    secant = [Fraction(1)]
    for degree in range(1, maximum + 1):
        secant.append(-sum(cosine[index] * secant[degree - index] for index in range(1, degree + 1)))
    secant5 = power(secant, 5)[:maximum + 1]
    polynomial = secant[:]
    secant_values = []
    for index in range(maximum + 1):
        value = polynomial[0]
        expected = general_recurrence(7, 1, index // 2) if index % 2 == 0 else 0
        assert value == expected
        secant_values.append(int(value))
        if index < maximum:
            polynomial = multiply(secant5, derivative(polynomial))[:maximum - index]
    factorial_cases = []
    for column in range(1, 5):
        polynomial = power(cosine if column < 3 else secant, abs(column - 3))[:maximum + 1]
        secant3 = power(secant, 3)[:maximum + 1]
        for index in range(maximum + 1):
            expected = general_recurrence(5, column, index // 2) if index % 2 == 0 else 0
            assert polynomial[0] == expected
            if index < maximum:
                polynomial = derivative(multiply(secant3, polynomial)[:maximum - index + 1])
        factorial_cases.append({'column': column, 'derivative_indices': [0, maximum], 'equal': True})
    return {'A281436_secant_formula': {'derivative_indices': [0, maximum], 'values': secant_values},
            'A281428_factorial_identity': factorial_cases,
            'normalization_witness': {'column': 1, 'index': 2, 'nested_derivative': 1, 'ordinary_coefficient': '1/2'},
            'offset_witnesses': [{'sequence': 'A281427', 'index': 1, 'sequence_value': 1, 'f_2_at_zero': 5},
                                 {'sequence': 'A281435', 'index': 1, 'sequence_value': 1, 'f_2_at_zero': 7}]}


def inverse_representation_checks():
    first = multiply([-1, 1], power([1, 2], 2))
    second = multiply([-1, 1], power([3, 4, 8], 2))
    assert first == [-1, -3, 0, 4]
    assert second == [-9, -15, -40, 0, 0, 64]
    maximum = 13
    sine = [Fraction(0)] * (maximum + 1)
    cosine = [Fraction(0)] * (maximum + 1)
    for index in range(7):
        sine[2 * index + 1] = Fraction(recurrence(2, 4, index), factorial(2 * index + 1))
        cosine[2 * index] = Fraction(recurrence(2, 1, index), factorial(2 * index))
    quotient = []
    for degree in range(maximum + 1):
        quotient.append(sine[degree] - sum(cosine[index] * quotient[degree - index]
                                         for index in range(1, degree + 1)))
    for index in range(7):
        assert quotient[2 * index + 1] == Fraction(recurrence(2, 1, index), factorial(2 * index + 1))
    return {'algebraic_eliminations': {'a_5_coefficients_in_C_squared': first,
                                      'a_7_coefficients_in_C_squared': second},
            'T_4_equals_integral_C_4_through_degree': maximum,
            'T_4_x5_coefficient': str(quotient[5]), 'printed_6_factorial_coefficient': '13/720'}


def second_column_period_checks():
    results = []
    for modulus in range(2, 16):
        remaining = 2 * modulus
        prime = 2
        local_bounds = []
        totient = remaining
        while prime * prime <= remaining:
            if remaining % prime == 0:
                totient = totient // prime * (prime - 1)
                height = 0
                while remaining % prime == 0:
                    height += 1
                    remaining //= prime
                local_bounds.append(prime ** (2 * height - 1))
            prime += 1
        if remaining > 1:
            totient = totient // remaining * (remaining - 1)
            local_bounds.append(remaining)
        onset = (max(local_bounds) + 2) // 2
        period = lcm(modulus, totient ** 2 // 2)
        values = column_residues(2, 4, 2 * modulus, onset + period + 10)
        observations = [(values[index] // 2 - 2 * index * values[index - 1]) % modulus
                        for index in range(onset, onset + period + 10)]
        assert all(observations[index + period] == observations[index] for index in range(10))
        results.append({'modulus': modulus, 'period': period, 'onset': onset,
                        'compared_indices': [onset, onset + 9], 'equal': True})
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    arguments = parser.parse_args()
    result = {'schema_version': 1, 'status': 'pass', 'arithmetic': 'exact integers and integer residues',
              'coefficient_checks': coefficient_checks(), 'certificates': exact_period_certificates(),
              'polynomial_operator_checks': polynomial_operator_checks(),
              'quadratic_operator_checks': quadratic_operator_checks(),
              'neighboring_coefficients': neighboring_coefficients(),
              'residue_row_certificates': residue_row_certificates(),
              'triangle_checks': triangle_checks(),
              'derivative_formula_checks': derivative_formula_checks(),
              'inverse_representation_checks': inverse_representation_checks(),
              'second_column_period_checks': second_column_period_checks()}
    content = json.dumps(result, indent=2, sort_keys=True) + '\n'
    if arguments.output:
        arguments.output.write_text(content)
    print(content, end='')


if __name__ == '__main__':
    main()
