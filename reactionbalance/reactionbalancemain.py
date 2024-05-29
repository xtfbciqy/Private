# -*- coding: utf-8 -*-
"""
Created on Tue May 28 12:53:05 2024

@author: 2acffb24
"""
from fractions import Fraction
import chemparse


def dimensions_func(matrix_para2):
    row_number_result = len(matrix_para2)
    col_number_result = len(matrix_para2[0])
    for i in range(1, row_number_result):
        if len(matrix_para2[i]) != col_number_result:
            return [row_number_result]
    return [row_number_result, col_number_result]


def identity_matrix_func(n_para2=-1, mn_para2=[-1]):
    if n_para2 >= 1:
        row_range = col_range = range(n_para2)
    elif mn_para2[0] >= 1:
        row_range = range(mn_para2[0])
        col_range = range(mn_para2[1])
    return [[1 if i == j else 0
             for j in col_range]
            for i in row_range]


def rationalize_func(matrix_para2, row_number_para2=-1, col_number_para2=-1):
    if row_number_para2 == -1:
        row_range = range(len(matrix_para2))
    else:
        row_range = range(row_number_para2)
    if col_number_para2 == -1:
        return [[Fraction(matrix_para2[i][j])
                 for j in range(len(matrix_para2[i]))]
                for i in row_range]
    else:
        col_range = range(col_number_para2)
        return [[Fraction(matrix_para2[i][j])
                 for j in col_range]
                for i in row_range]


def transpose_func(matrix_para2, row_number_para2=-1, col_number_para2=-1):
    if row_number_para2 == -1:
        row_range = range(len(matrix_para2))
    else:
        row_range = range(row_number_para2)
    if col_number_para2 == -1:
        col_range = range(len(matrix_para2[0]))
    else:
        col_range = range(col_number_para2)
    return [[matrix_para2[i][j]
             for i in row_range]
            for j in col_range]


def matrix_preprocessor_func(matrix_para1):
    dimensions = dimensions_func(matrix_para1)
    if len(dimensions) == 2:
        row_number = dimensions[0]
        col_number = dimensions[1]
        rc_sum = row_number + col_number
        return [transpose_func(
            rationalize_func(matrix_para1 + identity_matrix_func(col_number),
                             rc_sum, col_number), rc_sum, col_number),
                col_number, row_number]


def kernel_basis_func(p_output_para1):
    mat_p = p_output_para1[0]
    row_a = p_output_para1[1]
    col_a = p_output_para1[2]
    col_b = row_a + col_a
    piv_r = 1
    piv_c = 1
    while piv_r <= row_a and piv_c <= col_a:
        arg_m = piv_r
        for ite_r in range(piv_r, row_a + 1):
            if (abs(mat_p[ite_r - 1][piv_c - 1]) >
                    abs(mat_p[arg_m - 1][piv_c - 1])):
                arg_m = ite_r
        if mat_p[arg_m - 1][piv_c - 1] == 0:
            piv_c += 1
        else:
            mat_p[piv_r - 1], mat_p[arg_m - 1] =\
                mat_p[arg_m - 1], mat_p[piv_r - 1]
            for ite_r in range(piv_r + 1, row_a + 1):
                if mat_p[ite_r - 1][piv_c - 1] != 0:
                    f = (mat_p[ite_r - 1][piv_c - 1]
                         / mat_p[piv_r - 1][piv_c - 1])
                    mat_p[ite_r - 1][piv_c - 1] = Fraction(0)
                    for ite_c in range(piv_c + 1, col_b + 1):
                        mat_p[ite_r - 1][ite_c - 1] -= (
                            mat_p[piv_r - 1][ite_c - 1] * f)
            piv_r += 1
            piv_c += 1
    kernel_basis_result = []
    piv_r = 1
    piv_c = 1
    while piv_r <= row_a and piv_c <= col_a:
        for ite_c in range(piv_c, col_a + 1):
            if mat_p[piv_r - 1][ite_c - 1] != 0:
                piv_r += 1
        kernel_basis_result.append(mat_p[piv_r - 1][col_a:])
        piv_r += 1
    return kernel_basis_result


def solution_sum_func(basis_para2, n_para2, m_para2, coeffs_para2):
    sum_result = [Fraction(0)] * m_para2
    for i in range(n_para2):
        for j in range(m_para2):
            sum_result[j] += coeffs_para2[i] * basis_para2[i][j]
    return [sum_result, coeffs_para2, sum(sum_result)]


def generate_permutations_func(size_range_para2, n_para2):
    def recursive_func(current_para3):
        if len(current_para3) == n_para2:
            permutation_result.append(current_para3)
            return
        for i in size_range_para2:
            recursive_func(current_para3 + [i])
    permutation_result = []
    recursive_func([])
    return permutation_result


def list_multiplication_func(list1_para2, list2_para2, n_para2):
    mult_result = [0] * n_para2
    for i in range(n_para2):
        mult_result[i] += list1_para2[i] * list2_para2[i]
    return mult_result


def find_positive_integer_solution(basis_para1):
    m = len(basis_para1[0])
    n = len(basis_para1)
    if n == 1:
        best_solution = basis_para1[0].copy()
        while True:
            max_denominator = 0
            for i in range(m):
                if best_solution[i].denominator > max_denominator:
                    max_denominator = best_solution[i].denominator
            if max_denominator == 1:
                return [int(c) for c in best_solution]
            for i in range(m):
                best_solution[i] *= max_denominator
    max_sum = Fraction(0)
    for i in range(n):
        max_sum += abs(sum(basis_para1[i]))
    max_sum = int(max_sum * max_sum.denominator * 2)
    best_solution = []
    coeffs_previous_list = []
    while len(best_solution) == 0:
        min_sum = float('inf')
        coeffs_list = generate_permutations_func(range(max_sum + 1), n)
        if len(coeffs_previous_list) > 0:
            for coeffs in coeffs_previous_list:
                if coeffs in coeffs_list:
                    coeffs_list.remove(coeffs)
        coeffs_symbol_list = generate_permutations_func(range(-1, 2, 2), n)
        for coeffs in coeffs_list:
            for symbol in coeffs_symbol_list:
                coeffs_ws = list_multiplication_func(coeffs, symbol, n)
                candidate = solution_sum_func(basis_para1, n, m, coeffs_ws)
                if all(c > 0 and c.denominator == 1 for c in candidate[0]):
                    if candidate[2] < min_sum:
                        min_sum = candidate[2]
                        best_solution.clear()
                        best_solution.append(candidate[0])
                    elif candidate[2] == min_sum:
                        best_solution.append(candidate[0])
        if len(best_solution) > 0:
            return [int(c) for c in best_solution[0]]
        else:
            coeffs_previous_list = coeffs_list.copy()
            max_sum *= 2


def main(equation_para0):
    elements = []
    index = equation_para0.find(" → ")
    if index != -1:
        l_equation = equation_para0[:index]
        r_equation = equation_para0[index + len(" → "):]
    l_dict = [chemparse.parse_formula(component)
              for component
              in l_equation.replace(" + ", " ").split()]
    r_dict = [chemparse.parse_formula(component)
              for component
              in r_equation.replace(" + ", " ").split()]
    for iter_dict in r_dict:
        for element in iter_dict.keys():
            iter_dict[element] *= -1
    dict_c = l_dict + r_dict
    for iter_dict in dict_c:
        for element in iter_dict.keys():
            if element not in elements:
                elements.append(element)
    matrix = [[int(iter_dict.get(element, 0))
               for iter_dict in dict_c]
              for element in sorted(elements)]
    p_output = matrix_preprocessor_func(matrix)
    basis = kernel_basis_func(p_output)
    result = find_positive_integer_solution(basis)
    print("Result:")
    len_l = len(l_equation.replace(" + ", " ").split())
    len_r = len(r_equation.replace(" + ", " ").split())
    for iter_int in range(len(result)):
        if result[iter_int] == 1:
            print_result = ""
        else:
            print_result = str(result[iter_int]) + " "
        if iter_int < len_l - 1:
            print(f"{print_result}\
{l_equation.replace(" + ", " ").split()[iter_int]}", end=" + ")
        elif iter_int == len_l - 1:
            print(f"{print_result}\
{l_equation.replace(" + ", " ").split()[iter_int]}", end=" → ")
        elif iter_int < len_l + len_r - 1:
            print(f"{print_result}\
{r_equation.replace(" + ", " ").split()[iter_int - len_l]}",
             end=" + ")
        elif iter_int == len_l + len_r - 1:
            print(f"{print_result}\
{r_equation.replace(" + ", " ").split()[iter_int - len_l]}",
             end="")


if __name__ == '__main__':
    # equation = "MnS + As2Cr10O35 + H2SO4 → HMnO4 + AsH3 + CrS3O12 + H2O"
    # equation = "(Cr(N2H4CO)6)4(Cr(CN)6)3 + KMnO4 + H2SO4 → K2Cr2O7 + MnSO4 + CO2 + KNO3 + K2SO4 + H2O"
    # equation = "Cu+ + Fe → Fe+3 + Cu"
    equation = input("\nInput example: CO + CO2 + H2 → CH4 + H2O\n             : ")
    main(equation)
