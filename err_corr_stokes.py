import numpy as np
import pandas as pd
import cvxpy as cp
import openpyxl as op
from qutip import Bloch, Qobj
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys


def bloch_from_dataframe(df, axes):
    """
    Plot Bloch sphere
    :param df: pd.DataFrame
    :param axes:
    :return:
    """
    bloch = Bloch(fig=fig,axes=axes)

    bloch.vector_color = plt.get_cmap('viridis')(np.linspace(0, 1, len(df)))

    bloch.add_states([
        Qobj(
            [[row.Jxx, row.beta + 1j * row.gamma], [row.beta - 1j * row.gamma, row.Jyy]]
        ) for _, row in df.iterrows()
    ])

    bloch.make_sphere()

if __name__=="__main__":

    '''
	for i in range(1, len(sys.argv), 2):
		argument = sys.argv[i]
		argument_value = sys.argv[i + 1]
		if argument == '--filename':
			filename = argument_value
    '''
    #######################################################################
    #
    #   Want to minimize E_i = (I_i - Tr(L_i * J * L_i^\dag))^2
    #
    #   Subject to J >= 0 and Hermitian
    #
    #######################################################################

    #Define the optical element matrices
    l1 = np.array([[1, 0], [0, 0]])
    l2 = np.array([[0, 0], [0, 1]])
    l3 = 0.5 * np.array([[1, 1], [1, 1]])
    l4 = 0.5 * np.exp(-1j* np.pi / 4) * np.array([[1, 1j], [1, 1j]])
    l = [l1, l2, l3, l4]

    #Read in data from excel sheet and normalize intensities
    filename = 'HWP_polarimeter_low.xlsx'
    data = pd.read_excel(filename, sheet_name='measured')
    S = np.array([data.values[i][1:5] for i in range(data.values.shape[0])])
    for s in S:
        s /= (s[0])

    #Set up variables for optimization
    #   Each element of J is the coherency matrix for different polarization angles from the file
    J = [cp.Variable((2, 2), PSD=True) for i in range(data.values.shape[0])]

    #Apply constraints of positivity and hermitian to J
    cost = [(S[i][0] - (J[i][0][0] + J[i][1][1])) ** 2 + \
            (S[i][1] - (J[i][0][0] - J[i][1][1])) ** 2 + \
            (S[i][2] - 2 * np.real(J[i][0][1])) ** 2 + \
            (S[i][3] - 2 * np.imag(J[i][0][1])) ** 2
            for i in range(data.values.shape[0])]

    
    prob = [cp.Problem(cp.Minimize(cost[i])) for i in range(data.values.shape[0])]

    for i in range(data.values.shape[0]):
        prob[i].solve()
        #print(J[i].value)

    #Write the coherency matrix entries to the excel spreadsheet
	#First construct the dataframe
    columns = ['theta', 'Jxx', 'Jyy', 'beta', 'gamma', 'trace']

    temp = []
    for i in range(data.values.shape[0]):
	    temp.append([data.values[i][0], J[i].value[0][0], J[i].value[1][1], np.real(J[i].value[0][1]),\
                     np.imag(J[i].value[0][1]), np.trace(J[i].value.dot(J[i].value))])
    temp = np.array(temp)
    results = pd.DataFrame(temp, columns=columns, dtype=float)

    wb = pd.read_excel(filename, sheet_name=None)
    wb['Adjusted2'] = results
    with pd.ExcelWriter(filename) as writer:
    	for key in wb:
    		wb[key].to_excel(writer, sheet_name=key, index=False)

    ####################################################################################

    #Plot on Bloch sphere
    fig = plt.figure(figsize=plt.figaspect(0.3))
    fig_title = "Correction for " + filename
    fig.suptitle(fig_title)

    #First get results from measurement
    ax = fig.add_subplot(121, projection='3d', azim=0, elev=20)
    ax.set_title("Calculated from data", loc='left')

    bloch_from_dataframe(pd.read_excel(filename, sheet_name='calculated'), ax)

    #Now get adjusted results
    ax = fig.add_subplot(122, projection='3d', azim=0, elev=20)
    ax.set_title("Adjusted", loc='left')

    bloch_from_dataframe(pd.read_excel(filename, sheet_name='Adjusted2'), ax)

    #Set some plotting configs
    norm = mpl.colors.Normalize(vmin=0, vmax=90)
    sm = mpl.cm.ScalarMappable(cmap=plt.get_cmap('viridis'), norm=norm)
    sm.set_array([])
    cbaxes = fig.add_axes([0.9, 0.1, 0.03, 0.8])
    fig.colorbar(sm, ticks=np.linspace(0, 90, 10), cax=cbaxes)
             #boundaries=np.arange(-0.05,2.1,.1))

    plt.show()