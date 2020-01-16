import numpy as np
import pandas as pd
import cvxpy as cp
import matplotlib.pyplot as plt

def dop(filename, sheet_name):
    """
    Calculates the degree of polarization from data provided
    :param filename: the name of the data file (should be an excel file)
    :param sheet_name:  the name of the specific sheet
    """

    #Read in the data and slice out relevant elements
    data = pd.read_excel(filename, sheet_name=sheet_name)
    size = data.values.shape[0]
    rows = [data.values[i][1:5] for i in range(size)]
    theta = [data.values[i][0] for i in range(size)]

    #Calculate DOP for each row
    DOP = []
    for i in range(size):
        s0 = rows[i][0] + rows[i][1]
        s1 = rows[i][0] - rows[i][1]
        s2 = 2 * rows[i][2]
        s3 = 2 * rows[i][3]
        DOP.append(np.sqrt(s1**2 + s2**2 + s3**2)/s0)

    return theta, DOP

def plot_dop(filename, label, theta, dop, fig_num=1):
    """
    Function for plotting the dop vs theta
    :param theta: the degree the polarizer was set to for the trial
    :param dop:   the dop calculated by dop above
    :param fig_num: the number of the figure to avoid clashing definitions
    :param title: the title for the plot
    """
    fig = plt.figure(fig_num)
    ax = fig.add_subplot(111)
    ax.set_title(title)
    ax.set_xlabel('theta')
    ax.set_ylabel('DOP')

    ax.plot(theta, dop, label=label)
    
def max_dop(filename, eps):
    """
    Function to maximize DOP within an error tolerance
    """
    #Read in calculated coherency matrix elements and reconstruct matrices
    data = pd.read_excel(filename, sheet_name='calculated')
    size = data.values.shape[0]
    J_elems = [data.values[i][1:5] for i in range(size)]
    J = []
    for i in range(size):
        temp = np.array([[J_elems[i][0], J_elems[i][2] + 1j*J_elems[i][3]],
                         [J_elems[i][2] - 1j*J_elems[i][3], J_elems[i][1]]])
        J.append(temp)

    dop = []
    for i in range(size):
        rho = cp.Variable((2, 2), PSD=True)
        constraints = [cp.norm(J[i] - rho) <= eps, cp.trace(rho) <= 1]
        cost = cp.sqrt(2*cp.trace(rho) - 0.5)

        prob = cp.Problem(cp.Maximize(cost), constraints)
        prob.solve(qcp=True, solver='CVXOPT')
        
        s0 = rho.value[0][0] + rho.value[1][1]
        s1 = rho.value[0][0] - rho.value[1][1]
        s2 = 2 * np.real(rho.value[0][1])
        s3 = 2 * np.imag(rho.value[0][1])

        dop.append(np.sqrt(s1**2 + s2**2 + s3**2)/s0)

    theta = [data.values[i][0] for i in range(size)]
    fig = plt.figure(5)
    ax = fig.add_subplot(111)
    ax.set_title("Max DOP for " + filename)
    ax.set_xlabel("theta")
    ax.set_ylabel("DOP")
    ax.plot(theta, dop)

if __name__ == '__main__':
    filename_list = ['HWP_hand_high.xlsx', 'HWP_polarimeter_high.xlsx', 'HWP_hand_low.xlsx', 'HWP_polarimeter_low.xlsx']

    fig_num = 1
    for filename in filename_list:
        theta, calc_dop = dop(filename, 'calculated')
        rho_dop = dop(filename, 'rho')[1]

        fig = plt.figure(fig_num)
        ax = fig.add_subplot(111)
        ax.set_title(('DOP for ' + filename))
        ax.set_xlabel('theta')
        ax.set_ylabel('DOP')
        ax.plot(theta, calc_dop, label='calculated')
        ax.plot(theta, rho_dop, label='rho')
        ax.legend(numpoints=1)
        fig_num += 1

    #max_dop(filename, 1e-1)
    plt.show()