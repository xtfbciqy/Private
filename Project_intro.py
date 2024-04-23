# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 21:38:35 2024
version = 001a3
@author: 2acffb24
"""

import matplotlib
import chemparse

str__equation = input("\nInput example: CH4 + O2 → CO2 + H2O\n             \
: ")
# str__equation = "MnS + As2Cr10O35 + H2SO4 → \
# HMnO4 + AsH3 + CrS3O12 + H2O"
# str__equation = "CH4 + O2 → CO2 + H2O"

list__elements = []
int__index = str__equation.find(" → ")
if int__index != -1:
    str__l_equation = str__equation[:int__index]
    str__r_equation = str__equation[int__index + len(" → "):]
list__l_dict = [chemparse.parse_formula(str__component)
                for str__component
                in str__l_equation.replace(" + ", " ").split()]
list__r_dict = [chemparse.parse_formula(str__component)
                for str__component
                in str__r_equation.replace(" + ", " ").split()]
for dict__iter in list__r_dict:
    for str__element in dict__iter.keys():
        dict__iter[str__element] *= -1
list__dict = list__l_dict + list__r_dict
for dict__iter in list__dict:
    for str__element in dict__iter.keys():
        if str__element not in list__elements:
            list__elements.append(str__element)
list__output = [[str__element] + [int(dict__iter.get(str__element, 0))
                                  for dict__iter in list__dict]
                for str__element in sorted(list__elements)]

print("")
list__c_matrix = list__output
for list__line in list__c_matrix:
    for var__iter in list__line:
        if not isinstance(var__iter, int):
            print("{0:2}".format(var__iter), end=":(  ")
        else:
            print("{0:4g}".format(var__iter), end="   ")
    print("  )")
    del list__line[:1]
print("\n", end="")


def fun__gcd(int__a, int__b):
    if int__b == 0:
        return int__a
    return fun__gcd(int__b, int__a % int__b)


def fun__gcd_list(list__a):
    list__nzl = [abs(int__i) for int__i in list__a if int__i != 0]
    if len(list__nzl) == 0:
        return 1
    int__r = list__nzl[0]
    for int__t in list__nzl[1:]:
        int__r = fun__gcd(int__r, int__t)
    return int__r


def fun__lcm_list(list__a):
    list__nzl = [abs(int__i) for int__i in list__a if int__i != 0]
    if len(list__nzl) == 0:
        return 1
    int__r = list__nzl[0]
    for int__t in list__nzl[1:]:
        int__r = int__r * int__t // fun__gcd(int__r, int__t)
    return int__r


list__backup = list__c_matrix
int__rindex = len(list__c_matrix)
for list__iter in list__c_matrix:
    for int__iter in range(int__rindex):
        if len(list__iter) != len(list__c_matrix[int__iter]):
            raise ValueError("Error message.")
int__cindex = len(list__c_matrix[0])

int__pivot_r = 1
int__pivot_c = 1
while (int__pivot_r <= int__rindex and int__pivot_c <= int__cindex):
    int__argmax = int__pivot_r
    for int__r_iter in range(int__pivot_r, int__rindex + 1):
        if (abs(list__c_matrix[int__r_iter - 1][int__pivot_c - 1]) >
                abs(list__c_matrix[int__argmax - 1][int__pivot_c - 1])):
            int__argmax = int__r_iter
    if list__c_matrix[int__argmax - 1][int__pivot_c - 1] == 0:
        int__pivot_c += 1
    else:
        list__c_matrix[int__pivot_r - 1], list__c_matrix[int__argmax - 1] = \
            list__c_matrix[int__argmax - 1], list__c_matrix[int__pivot_r - 1]
        for int__r_iter in range(int__pivot_r + 1, int__rindex + 1):
            if list__c_matrix[int__r_iter - 1][int__pivot_c - 1] != 0:
                int__o1 = list__c_matrix[int__pivot_r - 1][int__pivot_c - 1]
                float__o2 = list__c_matrix[int__r_iter - 1][int__pivot_c - 1]\
                    / list__c_matrix[int__pivot_r - 1][int__pivot_c - 1]
                list__c_matrix[int__r_iter - 1][int__pivot_c - 1] = 0
                for int__c_iter in range(int__pivot_c + 1, int__cindex + 1):
                    list__c_matrix[int__r_iter - 1][int__c_iter - 1] = int(
                        round(int__o1 *
                              (list__c_matrix[int__r_iter - 1]
                               [int__c_iter - 1]
                               - float__o2 *
                               list__c_matrix[int__pivot_r - 1]
                               [int__c_iter - 1]), 0))
        for list__iter in list__c_matrix:
            int__gcm = fun__gcd_list(list__iter)
            if int__gcm > 1:
                for int__iter in range(int__cindex):
                    list__iter[int__iter] = list__iter[int__iter] // int__gcm
        int__pivot_r += 1
        int__pivot_c += 1

int__pivot_r = int__rindex
int__pivot_c = int__cindex - 1

while (int__pivot_r > 0 and int__pivot_c > 0):
    if list__c_matrix[int__pivot_r - 1][int__pivot_c - 1] == 0:
        int__pivot_r -= 1
    else:
        for int__r_iter in range(int__pivot_r - 1, 0, -1):
            if list__c_matrix[int__r_iter - 1][int__pivot_c - 1] != 0:
                int__o1 = list__c_matrix[int__pivot_r - 1][int__pivot_c - 1]
                float__o2 = list__c_matrix[int__r_iter - 1][int__pivot_c - 1]\
                    / list__c_matrix[int__pivot_r - 1][int__pivot_c - 1]
                for int__c_iter in range(1, int__cindex + 1):
                    list__c_matrix[int__r_iter - 1][int__c_iter - 1] = int(
                        round(
                            int__o1 *
                            (list__c_matrix[int__r_iter - 1]
                             [int__c_iter - 1] -
                             float__o2 *
                             list__c_matrix[int__pivot_r - 1]
                             [int__c_iter - 1]), 0))
        for list__iter in list__c_matrix:
            int__gcm = fun__gcd_list(list__iter)
            if int__gcm > 1:
                for int__iter in range(int__cindex):
                    list__iter[int__iter] = list__iter[int__iter] // int__gcm
        int__pivot_r -= 1
        int__pivot_c -= 1

list__lcm = []
for int__iter in range(1, int__rindex + 1):
    if int__iter <= int__cindex - 1:
        list__lcm.append(list__c_matrix[int__iter - 1][int__iter - 1])
int__lcm = fun__lcm_list(list__lcm)
for int__iter in range(1, int__rindex + 1):
    if (int__iter <= int__cindex - 1 and
        abs(list__c_matrix[int__iter - 1][int__iter - 1])
            != int__lcm and list__c_matrix[int__iter - 1][int__iter - 1] != 0):
        float__o3 = int__lcm / list__c_matrix[int__iter - 1][int__iter - 1]
        for int__c_iter in range(1, int__cindex + 1):
            list__c_matrix[int__iter - 1][int__c_iter - 1] = int(round(
                float__o3 * list__c_matrix[int__iter - 1][int__c_iter - 1], 0))
    if (int__iter <= int__cindex - 1 and
            list__c_matrix[int__iter - 1][int__iter - 1] > 0):
        for int__c_iter in range(1, int__cindex + 1):
            list__c_matrix[int__iter - 1][int__c_iter - 1] *= -1

list__result = []
for list__iter in list__c_matrix:
    print("      ", end="")
    for var__iter in list__iter:
        print("{0:4g}".format(var__iter), end="   ")
    print("\n", end="")
    if abs(list__iter[-1]) != 0:
        list__result.append(abs(list__iter[-1]))
list__result.append(abs(list__c_matrix[0][0]))
print("")
for list__iter in list__c_matrix:
    print("      ", end="")
    for var__iter in list__iter:
        if var__iter == 0:
            print("{:4g}".format(float(var__iter)), end="   ")
        elif var__iter != 0 and var__iter == list__c_matrix[0][0]:
            print("{:4g}".format(float(var__iter / list__c_matrix[0][0])),
                  end="   ")
        elif var__iter != list__c_matrix[0][0]:
            print("{}".format(float(var__iter / list__c_matrix[0][0])), end="")
    print("\n", end="")
print(f"\nKernel:\n\n{list__result}\n")
print("Result:\n")

int__len_l = len(str__l_equation.replace(" + ", " ").split())
int__len_r = len(str__r_equation.replace(" + ", " ").split())

for int__iter in range(len(list__result)):
    if list__result[int__iter] == 1:
        str__print_result = ""
    else:
        str__print_result = str(list__result[int__iter]) + " "
    if int__iter < int__len_l - 1:
        print(f"{str__print_result}\
{str__l_equation.replace(" + ", " ").split()[int__iter]}", end=" + ")
    elif int__iter == int__len_l - 1:
        print(f"{str__print_result}\
{str__l_equation.replace(" + ", " ").split()[int__iter]}", end=" → ")
    elif int__iter < int__len_l + int__len_r - 1:
        print(f"{str__print_result}\
{str__r_equation.replace(" + ", " ").split()[int__iter - int__len_l]}",
         end=" + ")
    elif int__iter == int__len_l + int__len_r - 1:
        print(f"{str__print_result}\
{str__r_equation.replace(" + ", " ").split()[int__iter - int__len_l]}",
         end="")
