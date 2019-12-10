import numpy as np
import pandas as pd
import cvxpy as cp

if __name__=="__main__":
	

	####################################################################################
	#
	#	Want to minimize the error expression:
	#		e(a1, a2, a3) = (I[i] - 0.5*Tr[l[i].(sigma0 + sigma[j]a[j]).l[i](\dag)])**2
	#
	#	This is a QCQP problem, so it needs to be cast as a SOCP
	#		min 	y + aTx
	#		s.t.	0.5 * ||x||**2 <= 1
	#				uTu <= y
	#				u = 1/sqrt(2) * LT.x
	#				0 <= y
	#		
	#		where
	#				LLT is the Cholesky decomposition of the parameter p0
	#
	#
	####################################################################################
	#Set up the vector variable to be optimized
	#This is a 3D column vector with components a1, a2, & a3
	x = cp.Variable(3)
	y = cp.Variable()
	u = cp.Variable(3)

	#Set up the parameters
	p0 = cp.Parameter((3,3), nonneg=True)
	p0.value = np.array([[4, 0, 0],
						 [0, 2, 0],
						 [0, 0, 2]])

	L = np.linalg.cholesky(p0.value)
	A = 1/np.sqrt(2) * np.identity(3)
	b = np.zeros(3)
	c = np.zeros(3)
	d = 1
	
	#Set up the constraints
	soc_constraints = [
      cp.SOC(c.T@x + d, A@x + b)
	]
	soc_constraints.append(cp.SOC(y, u))
	soc_constraints.append(y >= 0)

	#Read in data from excel sheet to construct q0
	#Just for theta = 0 right now
	data = pd.read_excel("11-19-19.xlsx", sheet_name='measured')
	I = data.values[0]
	I /= (I[0]+I[1])
	q0 = cp.Parameter(3)
	q0.value = np.array([I[1]-I[0], 0.5-I[2], 0.5-I[3]])

	#Create the objective and the problem
	objective = cp.Minimize(y + q0.T@x)
	prob = cp.Problem(objective, soc_constraints)
	
	
	prob.solve(solver=cp.CVXOPT)
	print("status:", prob.status)
	print("optimal x value", x.value)
	print("optimal y value", y.value)
	print("optimal u value", u.value)
	print("optimal var", x.value, y.value)
	
