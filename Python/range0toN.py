import numpy as np

def range0toN(A, myrange):
    a = myrange[0]
    b = myrange[1]
    # A[:] = [b if ele > b else a if ele < a else ele for ele in A]
    A[A > b] = b
    A[A < a] = a
    return A
