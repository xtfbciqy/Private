h := 1
k := 1

while h ≤ m and k ≤ n
    i_max := argmax (i = h ... m, abs(A[i, k]))
    if A[i_max, k] = 0
        k := k + 1
    else
        swap rows(h, i_max)
        for i = h + 1 ... m:
            f := A[i, k] / A[h, k]
            A[i, k] := 0
            for j = k + 1 ... n:
                A[i, j] := A[i, j] - A[h, j] * f
        h := h + 1
        k := k + 1

