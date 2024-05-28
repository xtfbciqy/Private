# -*- coding: utf-8 -*-
"""
Created on Sun May 26 16:43:20 2024

@author: 2acffb24
"""
from fractions import Fraction
import numpy as np
import matplotlib.pyplot as plt
import time
import reactionbalancemain


def generate_matrix_func(n=2, limit=10):
    np.random.seed(n)
    if n > 1:
        last_column = np.array(0)
        while 0 in last_column:
            kernel = np.array(0)
            while 0 in kernel:
                kernel = np.array([Fraction(np.random.uniform(0, limit))
                                   for _ in range(n)]).reshape(-1, 1)
            matrix_b = np.array(0)
            while 0 in matrix_b:
                matrix_b = np.array([[
                    Fraction(np.random.uniform(0, limit)) *
                    np.random.choice([-1, 1])
                    for _ in range(n - 1)] for _ in range(n)])
            last_column = -np.dot(matrix_b, kernel[:-1]) / kernel[-1]
        return np.hstack((matrix_b, last_column)).tolist()


def record_time(start_n=2, end_n=60, step_n=2):
    x = [0, 1]
    y = [0.0, 0.0]
    for n in range(start_n, end_n + 1, step_n):
        start_time = time.time()
        reactionbalancemain.find_positive_integer_solution(
            reactionbalancemain.kernel_basis_func(
                reactionbalancemain.matrix_preprocessor_func(
                    generate_matrix_func(n))))
        end_time = time.time() - start_time
        x.append(n)
        y.append(end_time)
    return x, y


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


def debug_func(n=60):
    func_input = generate_matrix_func(n)
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
    # debug_func()
    plot_figure(record_time())
