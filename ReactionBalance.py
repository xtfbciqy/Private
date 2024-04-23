# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 21:38:35 2024
version = 001a2
@author: 2acffb24
"""

import matplotlib
import chemparse

equation = input("\nInput example: CH4 + O2 → CO2 + H2O\n             \
: ")
# equation = "MnS + As2Cr10O35 + H2SO4 → \
# HMnO4 + AsH3 + CrS3O12 + H2O"
# equation = "CH4 + O2 → CO2 + H2O"

elements = []
index = equation.find(" → ")
if index != -1:
    l_equation = equation[:index]
    r_equation = equation[index + len(" → "):]
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
output = [[element] + [int(iter_dict.get(element, 0))
                       for iter_dict in dict_c]
          for element in sorted(elements)]

print("")
c_matrix = output
for line in c_matrix:
    for iter_var in line:
        if not isinstance(iter_var, int):
            print("{0:2}".format(iter_var), end=":(  ")
        else:
            print("{0:4g}".format(iter_var), end="   ")
    print("  )")
    del line[:1]
print("\n", end="")


def gcd(a, b):
    if b == 0:
        return a
    return gcd(b, a % b)


def gcd_list(a):
    nzl = [abs(i) for i in a if i != 0]
    if len(nzl) == 0:
        return 1
    r = nzl[0]
    for t in nzl[1:]:
        r = gcd(r, t)
    return r


def lcm_list(a):
    nzl = [abs(i) for i in a if i != 0]
    if len(nzl) == 0:
        return 1
    r = nzl[0]
    for t in nzl[1:]:
        r = r * t // gcd(r, t)
    return r


backup = c_matrix
rindex = len(c_matrix)
for iter_list in c_matrix:
    for iter_int in range(rindex):
        if len(iter_list) != len(c_matrix[iter_int]):
            raise ValueError("Error message.")
cindex = len(c_matrix[0])

pivot_r = 1
pivot_c = 1
while (pivot_r <= rindex and pivot_c <= cindex):
    argmax = pivot_r
    for r_iter in range(pivot_r, rindex + 1):
        if (abs(c_matrix[r_iter - 1][pivot_c - 1]) >
                abs(c_matrix[argmax - 1][pivot_c - 1])):
            argmax = r_iter
    if c_matrix[argmax - 1][pivot_c - 1] == 0:
        pivot_c += 1
    else:
        c_matrix[pivot_r - 1], c_matrix[argmax - 1] = \
            c_matrix[argmax - 1], c_matrix[pivot_r - 1]
        for r_iter in range(pivot_r + 1, rindex + 1):
            if c_matrix[r_iter - 1][pivot_c - 1] != 0:
                o1 = c_matrix[pivot_r - 1][pivot_c - 1]
                o2 = c_matrix[r_iter - 1][pivot_c - 1]\
                    / c_matrix[pivot_r - 1][pivot_c - 1]
                c_matrix[r_iter - 1][pivot_c - 1] = 0
                for c_iter in range(pivot_c + 1, cindex + 1):
                    c_matrix[r_iter - 1][c_iter - 1] = int(
                        round(o1 *
                              (c_matrix[r_iter - 1]
                               [c_iter - 1]
                               - o2 *
                               c_matrix[pivot_r - 1]
                               [c_iter - 1]), 0))
        for iter_list in c_matrix:
            gcm = gcd_list(iter_list)
            if gcm > 1:
                for iter_int in range(cindex):
                    iter_list[iter_int] = iter_list[iter_int] // gcm
        pivot_r += 1
        pivot_c += 1

pivot_r = rindex
pivot_c = cindex - 1

while (pivot_r > 0 and pivot_c > 0):
    if c_matrix[pivot_r - 1][pivot_c - 1] == 0:
        pivot_r -= 1
    else:
        for r_iter in range(pivot_r - 1, 0, -1):
            if c_matrix[r_iter - 1][pivot_c - 1] != 0:
                o1 = c_matrix[pivot_r - 1][pivot_c - 1]
                o2 = c_matrix[r_iter - 1][pivot_c - 1]\
                    / c_matrix[pivot_r - 1][pivot_c - 1]
                for c_iter in range(1, cindex + 1):
                    c_matrix[r_iter - 1][c_iter - 1] = int(
                        round(
                            o1 *
                            (c_matrix[r_iter - 1]
                             [c_iter - 1] -
                             o2 *
                             c_matrix[pivot_r - 1]
                             [c_iter - 1]), 0))
        for iter_list in c_matrix:
            gcm = gcd_list(iter_list)
            if gcm > 1:
                for iter_int in range(cindex):
                    iter_list[iter_int] = iter_list[iter_int] // gcm
        pivot_r -= 1
        pivot_c -= 1

lcm = []
for iter_int in range(1, rindex):
    if iter_int <= cindex - 1:
        lcm.append(c_matrix[iter_int - 1][iter_int - 1])
lcm = lcm_list(lcm)
for iter_int in range(1, rindex):
    if (iter_int <= cindex - 1 and
        abs(c_matrix[iter_int - 1][iter_int - 1])
            != lcm and c_matrix[iter_int - 1][iter_int - 1] != 0):
        o3 = lcm / c_matrix[iter_int - 1][iter_int - 1]
        for c_iter in range(1, cindex + 1):
            c_matrix[iter_int - 1][c_iter - 1] = int(round(
                o3 * c_matrix[iter_int - 1][c_iter - 1], 0))
    if (iter_int <= cindex - 1 and
            c_matrix[iter_int - 1][iter_int - 1] > 0):
        for c_iter in range(1, cindex + 1):
            c_matrix[iter_int - 1][c_iter - 1] *= -1

result = []
for iter_list in c_matrix:
    print("      ", end="")
    for iter_var in iter_list:
        print("{0:4g}".format(iter_var), end="   ")
    print("\n", end="")
    if abs(iter_list[-1]) != 0:
        result.append(abs(iter_list[-1]))
result.append(abs(c_matrix[0][0]))
print("")
for iter_list in c_matrix:
    print("      ", end="")
    for iter_var in iter_list:
        if iter_var == 0:
            print("{:4g}".format(float(iter_var)), end="   ")
        if iter_var != 0:
            print("{:4g}".format(float(iter_var / c_matrix[0][0])), end="   ")
    print("\n", end="")
print(f"\nKernel:\n\n{result}\n")
print("Result:\n")

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
