# -*- coding: utf-8 -*-
from fractions import Fraction
import matplotlib.pyplot as plt
import random
import time
import reactionbalancemain


def generate_matrix(n=2, limit=10):
    random.seed(n)

    def uniform_fraction(limit):
        return Fraction(random.uniform(0, limit))
    if n > 1:
        last_column = [0]
        while 0 in last_column:
            kernel = [0]
            while 0 in kernel:
                kernel = [uniform_fraction(limit) for _ in range(n)]
            matrix_b = [[0]]
            while any(0 in row for row in matrix_b):
                matrix_b = [[uniform_fraction(limit) * random.choice([-1, 1])
                             for _ in range(n - 1)] for _ in range(n)]
            last_column = [-sum(matrix_b[i][j] * kernel[j]
                                for j in range(n - 1)) / kernel[-1]
                           for i in range(n)]
        for i in range(n):
            matrix_b[i].append(last_column[i])
        return matrix_b


def record_time(start_n=2, end_n=60, step_n=2):
    x = [0, 1]
    y = [0.0, 0.0]
    for n in range(start_n, end_n + 1, step_n):
        func_input = generate_matrix(n)
        start_time = time.time()
        reactionbalancemain.find_positive_integer_solution(
            reactionbalancemain.kernel_basis_func(
                reactionbalancemain.matrix_preprocessor_func(
                    func_input)))
        end_time = time.time() - start_time
        x.append(n)
        y.append(end_time)
    return [x, y]


def plot_figure(list_xy):
    plt.figure(figsize=(8, 5), dpi=480)
    plt.plot(list_xy[0], list_xy[1], color='blue', linestyle='-',
             linewidth=1.5, label='Line 1')
    plt.xlabel('Input Size n (n x n Matrix)')
    plt.ylabel('Operating Time (Second)')
    plt.grid(True)
    plt.legend()
    plt.savefig('figure.svg')
    plt.savefig('figure.png')


def debug(n=60):
    func_input = generate_matrix(n)
    result = reactionbalancemain.find_positive_integer_solution(
        reactionbalancemain.kernel_basis_func(
            reactionbalancemain.matrix_preprocessor_func(
                func_input)))
    print("{", end="")
    for i in range(len(func_input)):
        print("{", end="")
        for j in range(len(func_input[i])):
            print(f"\
({func_input[i][j].numerator}/{func_input[i][j].denominator})", end="")
            if j < len(func_input[i]) - 1:
                print(", ", end="")
        print("}", end="")
        if i < len(func_input) - 1:
            print(", ", end="")
        else:
            print("}", end="")
    print("\n\n{", end="")
    for i in range(len(result)):
        print(f"{result[i]}", end="")
        if i < len(func_input) - 1:
            print(", ", end="")
        else:
            print("}", end="")


if __name__ == '__main__':
    # debug()
    plot_figure(record_time())
