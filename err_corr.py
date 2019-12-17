import numpy as np
import pandas as pd
import cvxpy as cp
import openpyxl as op
from qutip import Bloch, Qobj
import matplotlib.pyplot as plt

fig = plt.figure(figsize=plt.figaspect(0.3))

def bloch_from_dataframe(df, axes):
    """
    Plot Bloch sphere
    :param df: pd.DataFrame
    :param axes:
    :return:
    """
    bloch = Bloch(axes=axes)

    bloch.vector_color = plt.cm.viridis(np.linspace(0, 1, len(df)))

    bloch.add_states([
        Qobj(
            [[row.Jxx, row.beta + 1j * row.gamma], [row.beta - 1j * row.gamma, row.Jyy]]
        ) for _, row in df.iterrows()
    ])

    bloch.make_sphere()

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
	soc_constraints.append(cp.norm(x) <= 1)

	#Read in data from excel sheet to construct q0
	#Just for theta = 0 right now
	filename = "11-19-19.xlsx"
	data = pd.read_excel(filename, sheet_name='measured')
	I = np.array([data.values[i][1:5] for i in range(data.values.shape[0])])
	q0_arr = []

	for i in I:
		#Normalize the intensities
		i /= (i[0] + i[1])/2
		q0 = cp.Parameter(3)
		q0.value = np.array([i[1]-i[0], 0.5-i[2], 0.5-i[3]])
		q0_arr.append(q0)

	q0_arr = np.array(q0_arr)

	#Create the objectives and the problems
	objectives = [cp.Minimize(y+q0_arr[i].T@x) for i in range(data.values.shape[0])]
	problems = []
	x_vals = []
	y_vals = []
	u_vals = []

	for i in range(len(objectives)):
		problems.append(cp.Problem(objectives[i], soc_constraints))
		problems[i].solve()
		x_vals.append(x.value)
		y_vals.append(y.value)
		u_vals.append(u.value)
		
	#######################################################################################

	#Check that the solution work	
	#	Initialize the sigmas
	sigma0 = np.identity(2)
	sigma1 = np.array([[1, 0],
					   [0, -1]])
	sigma2 = np.array([[0, 1],
					   [1, 0]])
	sigma3 = np.array([[0, 1j],
					   [-1j, 0]])
	sigmas = [sigma0, sigma1, sigma2, sigma3]
	
	#	Initialize the a vector
	a = []
	for i in range(len(x_vals)):
		a0 = 1 + 0j
		a1 = x_vals[i][0] + 0j
		a2 = x_vals[i][1] + 0j
		a3 = x_vals[i][2] + 0j
		a.append([a0, a1, a2, a3])

	J = np.zeros(shape=(data.values.shape[0], 2, 2), dtype=complex)
	for i in range(data.values.shape[0]):
		for j in range(len(sigmas)):
			J[i] += 0.5 * (a[i][j]*sigmas[j])

	#######################################################################################

	#Write the coherency matrix entries to the excel spreadsheet
	#First construct the dataframe
	columns = ['theta', 'Jxx', 'Jyy', 'beta', 'gamma', 'trace']
	
	temp = []
	for i in range(data.values.shape[0]):
		temp.append([data.values[i][0], J[i][0][0], J[i][1][1], np.real(J[i][0][1]), np.imag(J[i][0][1]), np.trace(J[i].dot(J[i]))])
	temp = np.array(temp)
	print(temp[0])
	results = pd.DataFrame(temp, columns=columns, dtype=float)

	wb = pd.read_excel(filename, sheet_name=None)
	wb['Adjusted'] = results
	with pd.ExcelWriter(filename) as writer:
		for key in wb:
			wb[key].to_excel(writer, sheet_name=key, index=False)
	
	

	#######################################################################################
	
	#Plot on Bloch sphere
	#First get results from measurement
	ax = fig.add_subplot(121, projection='3d', azim=0, elev=20)
	ax.set_title("Calculated from data")

	bloch_from_dataframe(pd.read_excel(filename, sheet_name='calculated'), ax)

	#Now get adjusted results
	ax = fig.add_subplot(122, projection='3d', azim=0, elev=20)
	ax.set_title("Adjusted")

	bloch_from_dataframe(pd.read_excel(filename, sheet_name='Adjusted'), ax)

	#Set some plotting configs
	import matplotlib as mpl
	norm = mpl.colors.Normalize(vmin=0, vmax=90)
	sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=norm)
	sm.set_array([])
	fig.colorbar(sm, ticks=np.linspace(0, 90, 10))
             #boundaries=np.arange(-0.05,2.1,.1))

	plt.show()