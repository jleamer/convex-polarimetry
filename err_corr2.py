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

def int_corr(filename):
    """
    Function to minimize the coherency matrix from the provided data
    :param data: pd.DataFrame object
    :return:     the list of minimized coherency matrices
    """
    wb = pd.read_excel(filename, sheet_name=None)
    data = wb['measured']
    size = data.values.shape[0]
    min_J = []
    columns = ['theta', 'Jxx', 'Jyy', 'beta', 'gamma', 'trace']

    #Slice the intensities out of the data file and normalize them
    I = np.array([data.values[i][1:5] for i in range(size)])
    for i in I:
        i /= (i[0] + i[1])

    for i in range(size):
        J = cp.Variable((2, 2), PSD=True)
        cost = cp.abs(I[i][0] -cp.trace((l1.conj().T).dot(l1) @ J)) ** 2 + \
               cp.abs(I[i][1] -cp.trace((l2.conj().T).dot(l2) @ J)) ** 2 + \
               cp.abs(I[i][2] -cp.trace((l3.conj().T).dot(l3) @ J)) ** 2 + \
               cp.abs(I[i][3] -cp.trace((l4.conj().T).dot(l4) @ J)) ** 2

        constraint = [cp.trace(J) <= 1]
        prob = cp.Problem(cp.Minimize(cost), constraint)
        prob.solve()
        
        min_J.append([data.values[i][0], J.value[0][0], J.value[1][1], np.real(J.value[0][1]), np.imag(J.value[0][1]), np.trace(J.value.dot(J.value))])

    min_J = np.array(min_J)
    temp = pd.DataFrame(min_J, columns=columns)
    wb['adjusted'] = temp

    with pd.ExcelWriter(filename) as writer:
        for key in wb:
            wb[key].to_excel(writer, sheet_name=key, index=False)

def stokes_corr(filename):
    """
    Function to minimize the coherency matrix from the provided data
    :param data: pd.DataFrame object
    :return:     the list of minimized coherency matrices
    """
    wb = pd.read_excel(filename, sheet_name=None)
    data = wb['measured']
    size = data.values.shape[0]
    min_J = []
    columns = ['theta', 'Jxx', 'Jyy', 'beta', 'gamma', 'trace']

    #Slice the Stokes parameters out of the data file and normalize them
    S = np.array([data.values[i][1:5] for i in range(size)])
    for s in S:
        s /= s[0] 

    for i in range(size):
        J = cp.Variable((2, 2), PSD=True)
        cost = (S[i][0] - (J[0][0] + J[1][1])) ** 2 + \
               (S[i][1] - (J[0][0] - J[1][1])) ** 2 + \
               (S[i][2] - 2 * np.real(J[0][1])) ** 2 + \
               (S[i][3] - 2 * np.imag(J[0][1])) ** 2

        constraint = [cp.trace(J) <= 1]
        prob = cp.Problem(cp.Minimize(cost), constraint)
        prob.solve()
        
        min_J.append([data.values[i][0], J.value[0][0], J.value[1][1], np.real(J.value[0][1]), np.imag(J.value[0][1]), np.trace(J.value.dot(J.value))])

    min_J = np.array(min_J)
    temp = pd.DataFrame(min_J, columns=columns)
    wb['adjusted'] = temp

    with pd.ExcelWriter(filename) as writer:
        for key in wb:
            wb[key].to_excel(writer, sheet_name=key, index=False)

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
    filename = 'HWP_hand_low.xlsx'
    int_corr(filename)

    filename2 = 'HWP_polarimeter_low.xlsx'
    stokes_corr(filename)
    

    ####################################################################################

    #Plot on Bloch sphere
    fig = plt.figure(1, figsize=plt.figaspect(0.3))
    fig_title = "Correction for " + filename
    fig.suptitle(fig_title)

    #First get results from measurement
    ax = fig.add_subplot(121, projection='3d', azim=0, elev=20)
    ax.set_title("Calculated from data", loc='left')

    bloch_from_dataframe(pd.read_excel(filename, sheet_name='calculated'), ax)

    #Now get adjusted results
    ax = fig.add_subplot(122, projection='3d', azim=0, elev=20)
    ax.set_title("Adjusted", loc='left')

    bloch_from_dataframe(pd.read_excel(filename, sheet_name='adjusted'), ax)

    #Set some plotting configs
    norm = mpl.colors.Normalize(vmin=0, vmax=90)
    sm = mpl.cm.ScalarMappable(cmap=plt.get_cmap('viridis'), norm=norm)
    sm.set_array([])
    cbaxes = fig.add_axes([0.9, 0.1, 0.03, 0.8])
    fig.colorbar(sm, ticks=np.linspace(0, 90, 10), cax=cbaxes)
             #boundaries=np.arange(-0.05,2.1,.1))

    fig2 = plt.figure(2, figsize=plt.figaspect(0.3))
    fig2_title = "Correction for " + filename2
    fig2.suptitle(fig2_title)

    ax2 = fig2.add_subplot(121, projection='3d', azim=0, elev=20)
    ax2.set_title("Calculated from data", loc='left')

    bloch_from_dataframe(pd.read_excel(filename2, sheet_name='calculated'), ax2)
    
    ax2 = fig2.add_subplot(122, projection='3d', azim=0, elev=20)
    ax2.set_title("Adjusted", loc='left')

    bloch_from_dataframe(pd.read_excel(filename2, sheet_name='adjusted'), ax2)
    cbaxes = fig2.add_axes([0.9, 0.1, 0.03, 0.8])
    fig2.colorbar(sm, ticks=np.linspace(0, 90, 10), cax=cbaxes)
    plt.show()