import cvxpy as cp
import numpy as np 

epsilon = 1e-1

rho = cp.Variable((2, 2), PSD=True)
J = np.random.rand(2,2)

cost = cp.norm(rho - J) ** 2

count = 0
while count < 10:
    try:
        constraints = [cp.trace(rho) == 1, cp.norm(rho - J) <= epsilon]

        prob = cp.Problem(cp.Minimize(cost), constraints)
        prob.solve()

        print(rho.value)
        break

    except:
        epsilon = epsilon/2