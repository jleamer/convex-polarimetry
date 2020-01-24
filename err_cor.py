import numpy as np
import pandas as pd
import cvxpy as cp
import openpyxl as op
from qutip import Bloch, Qobj
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys


def bloch_from_dataframe(df, fig, axes):
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

def plot_Bloch(filename, fig_num):
    """
    Function for plotting the Bloch spheres created by bloch_from_dataframe
    :param filename: name of the data file to be used for plotting (should be an excel file)
    :param fig_num:  the number to be assigned to the figure
    """
    #Plot on Bloch sphere
    fig = plt.figure(fig_num, figsize=plt.figaspect(0.3))
    fig_title = "Correction for " + filename
    fig.suptitle(fig_title)

    #First get results from measurement
    ax = fig.add_subplot(121, projection='3d', azim=0, elev=20)
    ax.set_title("Calculated from data", loc='left')

    bloch_from_dataframe(pd.read_excel(filename, sheet_name='calculated'), fig, ax)

    #Now get adjusted results
    ax = fig.add_subplot(122, projection='3d', azim=0, elev=20)
    ax.set_title("Adjusted", loc='left')

    bloch_from_dataframe(pd.read_excel(filename, sheet_name='rho'), fig, ax)

    #Set some plotting configs
    norm = mpl.colors.Normalize(vmin=0, vmax=90)
    sm = mpl.cm.ScalarMappable(cmap=plt.get_cmap('viridis'), norm=norm)
    sm.set_array([])
    cbaxes = fig.add_axes([0.9, 0.1, 0.03, 0.8])
    fig.colorbar(sm, ticks=np.linspace(0, 90, 10), cax=cbaxes)
             #boundaries=np.arange(-0.05,2.1,.1))

def compute_rho(filename):
    """
    Least squares minimize the calculated coherency matrix
    :param filename: the name of the data file (should be an excel file)
    """

    #Construct the dataframe from the file
    wb = pd.read_excel(filename, sheet_name=None)
    data = wb['calculated']
    size = data.values.shape[0]
    columns = ['theta', 'Jxx', 'Jyy', 'beta', 'gamma', 'trace_sq']

    #Slice coherency matrix elements from the data
    J_elems = [data.values[i][1:5] for i in range(size)]

    #Set up list for storing rhos
    rho_list = []
    temp = []

    for i in range(size):
        #Use J_elems to construct matrix
        J = np.array([[J_elems[i][0], J_elems[i][2] + 1j*J_elems[i][3]],
                      [J_elems[i][2] - 1j*J_elems[i][3], J_elems[i][1]]])

        #Construct variable to be minimized
        rho = cp.Variable((2, 2), PSD=True)

        #Construct constraint(s)
        constraints = [cp.trace(rho) == 1]

        #Construct cost function
        cost = cp.norm(J - rho) ** 2

        #Construct problem
        prob = cp.Problem(cp.Minimize(cost), constraints)
        prob.solve(solver='CVXOPT')

        rho_list.append(rho.value[0])

        #Construct a temporary list to be added to the excel file
        temp.append([data.values[i][0], rho.value[0][0], rho.value[1][1], np.real(rho.value[0][1]), np.imag(rho.value[0][1]), np.sum(rho.value.T*rho.value)])

    #Construct the new dataframe and write to the file
    results = pd.DataFrame(np.array(temp), columns=columns)
    wb['rho'] = results
    with pd.ExcelWriter(filename) as writer:
        for key in wb:
            wb[key].to_excel(writer, sheet_name=key, index=False)

def plot_error(filename):
    """
    """
    calc = pd.read_excel(filename, sheet_name='calculated')
    rho = pd.read_excel(filename, sheet_name='rho')
    size = calc.values.shape[0]

    J_elems = [calc.values[i][1:5] for i in range(size)]
    rho_elems = [rho.values[i][1:5] for i in range(size)]
    theta = [calc.values[i][0] for i in range(size)]

    error = []
    for i in range(size):
        J = np.array([[J_elems[i][0], J_elems[i][2] + 1j*J_elems[i][3]],
                      [J_elems[i][2] - 1j*J_elems[i][3], J_elems[i][1]]])
        rho = np.array([[rho_elems[i][0], rho_elems[i][2] + 1j*rho_elems[i][3]],
                        [rho_elems[i][2] - 1j*rho_elems[i][3], rho_elems[i][1]]])
        error.append(np.linalg.norm(J - rho))

    fig = plt.figure(max(plt.get_fignums()) + 1)
    ax = fig.add_subplot(111)
    ax.set_title('Error for ' + filename)
    ax.set_xlabel('theta')
    ax.set_ylabel('Error')
    ax.plot(theta, error)

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
    #   Want to minimize || J_calc - rho || ** 2
    #
    #   Subject to J >= 0 and trace(J) <= 1
    #
    #######################################################################

    #Define the optical element matrices
    l1 = np.array([[1, 0], [0, 0]])
    l2 = np.array([[0, 0], [0, 1]])
    l3 = 0.5 * np.array([[1, 1], [1, 1]])
    l4 = 0.5 * np.exp(-1j* np.pi / 4) * np.array([[1, 1j], [1, 1j]])
    l = [l1, l2, l3, l4]

    #Read in data from excel sheet and normalize intensities
    fig_num = 1
    filename_list = ['HWP_hand_high.xlsx', 'HWP_polarimeter_high.xlsx', 'HWP_hand_low.xlsx', 'HWP_polarimeter_low.xlsx']
    for filename in filename_list:
        compute_rho(filename)
        plot_Bloch(filename, fig_num)
        fig_num += 1

    for filename in filename_list:
        plot_error(filename)

    plt.show()