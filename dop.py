import numpy as np
import pandas as pd
import cvxpy as cp
from scipy.linalg import sqrtm
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
    dop = []
    for i in range(size):
        s0 = rows[i][0] + rows[i][1]
        s1 = rows[i][0] - rows[i][1]
        s2 = 2 * rows[i][2]
        s3 = 2 * rows[i][3]
        dop.append(np.sqrt(s1**2 + s2**2 + s3**2)/s0)

    return theta, dop


def adaptive_step(J, eps_guess):
    """
    """
    prev_eps = 0
    max_eps = 0
    epsilons = [eps_guess]
    eps = eps_guess
    count = 0

    rho = cp.Variable((2,2), PSD=True)
    cost = cp.norm(J - rho) ** 2

    #First find the largest value for eps that breaks CVXPY
    while True:
        try:
            constraints = [cp.trace(rho) == 1, cp.norm(J - rho) <= eps]
            prob = cp.Problem(cp.Minimize(cost), constraints)
            prob.solve(solver='CVXOPT')
            
            if rho.value is None:
                raise ValueError

            eps *= 2
        
        except:
            max_eps = eps
            eps = eps_guess
            break

    #Now perform bisection method to find the closest eps that breaks
    eps_mid = np.abs(max_eps - eps_guess)/2
    while count < 10:
        try:
            constraints = [cp.trace(rho) == 1, cp.norm(J - rho) <= eps_mid]
            prob = cp.Problem(cp.Minimize(cost), constraints)
            prob.solve(solver='CVXOPT')
            
            if rho.value is None:
                raise ValueError
            
            #if the code gets here, it succeeded and we can go higher
            eps_guess = eps_mid
            count += 1
            
        
        except:
            #if eps_mid fails, we need to go lower
            max_eps = eps_mid
            eps_mid = np.abs(max_eps - eps_guess)/2
            count += 1

    return rho, max_eps


def min_dop(filename, eps_guess, fig_num):
    """
    Function to maximize DOP within an error tolerance
    Each rho has its own optimal epsilon in this function -- maybe a problem?
    :param filename: the name of the file to maximize dop on
    :param eps: the initial error bound
    :param fig_num: the number of the figure to be plotted on
    """
    #Read in calculated coherency matrix elements and reconstruct matrices
    data = pd.read_excel(filename, sheet_name='calculated')
    size = data.values.shape[0]
    J_elems = [data.values[i][1:5] for i in range(size)]

    dop = []
    epsilons = []
    
    for i in range(size):
        prev_eps = 0
        J = np.array([[J_elems[i][0], J_elems[i][2] + 1j*J_elems[i][3]],
                      [J_elems[i][2] - 1j*J_elems[i][3], J_elems[i][1]]])
        
        rho = adaptive_step(J, eps_guess)[0]

        temp = rho.value
        s0 = temp[0][0] + temp[1][1]
        s1 = temp[0][0] - temp[1][1]
        s2 = 2 * np.real(temp[0][1])
        s3 = 2 * np.imag(temp[0][1])

        dop.append(np.sqrt(s1**2 + s2**2 + s3**2)/s0)

    theta = [data.values[i][0] for i in range(size)]
    fig = plt.figure(fig_num)
    ax = fig.add_subplot(111)
    ax.set_title("Max DOP for " + filename)
    ax.set_xlabel("theta")
    ax.set_ylabel("DOP")
    ax.plot(theta, dop, '--', label='Worst')
    ax.legend(numpoints=1)


def min_dop_same(filename, eps_guess, fig_num):
    """
    Function to maximize DOP within an error tolerance
    Each rho has its own optimal epsilon in this function -- maybe a problem?
    :param filename: the name of the file to maximize dop on
    :param eps: the initial error bound
    :param fig_num: the number of the figure to be plotted on
    """
    #Read in calculated coherency matrix elements and reconstruct matrices
    data = pd.read_excel(filename, sheet_name='calculated')
    size = data.values.shape[0]
    J_elems = [data.values[i][1:5] for i in range(size)]

    dop = []
    epsilons = []

    for i in range(size):
        J = np.array([[J_elems[i][0], J_elems[i][2] + 1j*J_elems[i][3]],
                      [J_elems[i][2] - 1j*J_elems[i][3], J_elems[i][1]]])

        epsilons.append(adaptive_step(J, eps_guess)[1])

    #opt_eps = max(epsilons)
    opt_eps = 3.2
    print(epsilons)
    print(opt_eps)

    for i in range(size):
        print(i)
        J = np.array([[J_elems[i][0], J_elems[i][2] + 1j*J_elems[i][3]],
                      [J_elems[i][2] - 1j*J_elems[i][3], J_elems[i][1]]])

        rho = cp.Variable((2, 2), PSD=True)
        cost = cp.norm(J - rho) ** 2
        constraints = [cp.trace(rho) == 1, cp.norm(J - rho) <= opt_eps]
        prob = cp.Problem(cp.Minimize(cost), constraints)
        prob.solve(solver='CVXOPT')

        temp = rho.value
        s0 = temp[0][0] + temp[1][1]
        s1 = temp[0][0] - temp[1][1]
        s2 = 2 * np.real(temp[0][1])
        s3 = 2 * np.imag(temp[0][1])

        dop.append(np.sqrt(s1**2 + s2**2 + s3**2)/s0)

    theta = [data.values[i][0] for i in range(size)]
    fig = plt.figure(fig_num)
    ax = fig.add_subplot(111)
    ax.set_title("Max DOP for " + filename)
    ax.set_xlabel("theta")
    ax.set_ylabel("DOP")
    ax.plot(theta, dop, label='Worst Same')


if __name__ == '__main__':
    filename_list = ['HWP_hand_high.xlsx', 'HWP_polarimeter_high.xlsx', 'HWP_hand_low.xlsx', 'HWP_polarimeter_low.xlsx']
    fig_num = 1
    '''
    fig = plt.figure(1)
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    fig2 = plt.figure(2)
    ax3 = fig2.add_subplot(121)
    ax4 = fig2.add_subplot(122)
    ax = [ax1, ax2, ax3, ax4]
    title = ['(a)', '(b)', '(c)', '(d)']

    ax[0].set_ylabel('DOP', fontsize=14)
    ax[2].set_ylabel('DOP', fontsize=14)
    for i in range(len(filename_list)):
        theta, calc_dop = dop(filename_list[i], 'calculated')
        rho_dop = dop(filename_list[i], 'rho')[1]

        ax[i].set_title(title[i])
        ax[i].set_xlabel(r'$\theta$' + ' (degrees)', fontsize=14)
        ax[i].plot(theta, calc_dop, label=r'$\mathbf{J}_{input}$', linewidth=2)
        ax[i].plot(theta, rho_dop, label=r'$\mathbf{J}_{unknown}$', linewidth=2)
        ax[i].legend(numpoints=1)
    '''
    
    titles = ['(a)', '(b)', '(c)', '(d)']
    for i in range(len(filename_list)):
        theta, calc_dop = dop(filename_list[i], 'calculated')
        rho_dop = dop(filename_list[i], 'rho')[1]

        fig = plt.figure(fig_num)
        ax = fig.add_subplot(111)
        ax.set_title(titles[i], fontsize=16)
        ax.set_xlabel(r'$\theta$' + ' (degrees)', fontsize=14)
        ax.set_ylabel('DOP', fontsize=14)
        ax.plot(theta, calc_dop, label=r'$\mathbf{J}_{input}$', linewidth=2)
        ax.plot(theta, rho_dop, label= r'$\mathbf{J}_{unknown}$', linewidth=2)
        ax.legend(numpoints=1)

        min_dop(filename_list[i], 1e-1, fig_num)
        fig_num += 1
    
    '''
    # Note that we found that rho isn't actually able to give back complex valued answers yet... Need to work on that
    # TODO: Ask question about performing optimization on hermitian variables and look into matlab version
    j0 = 0.5* np.array([[1, -1j], [1j, 1]])
    print(np.linalg.eigvalsh(j0))
    test = cp.Variable((2,2), hermitian=True)
    cost = cp.norm(test - j0)
    constraints = [cp.trace(test) <= 1]
    
    prob = cp.Problem(cp.Minimize(cost))
    prob.solve(solver='CVXOPT')
    print(test.value)
    print(np.linalg.eigvalsh(test.value))
    '''
    # TODO: Iterate over epsilon using bisection to find optimal epsilon value (worst case scenario)
    #min_dop(filename_list[0], 1e-1, 1)
    #min_dop_same(filename_list[0], 1e-1, 1)
    plt.show()
