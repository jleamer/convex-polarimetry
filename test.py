import cvxpy as cp
import numpy as np

A = np.random.normal(size=(2,2))
B = np.random.normal(size=(2,2))

print(np.trace(A.dot(B)))
print(np.sum(B.T * A))