from partial_rejection import *
from square_grid import create_grid_x_list, create_grid_bad_event_dict
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def grid_results_graph(m, trials, low, high, skip=1, p=0.5) :
    side_lengths = []
    iterations = []
    for i in range(low, high, skip):
        avg_iterations = 0
        side_lengths.append(i)
        x_list = create_grid_x_list(i, m, p)
        bad_event_dict = create_grid_bad_event_dict(i, m)
        for j in range(trials):
            avg_iterations += general_partial_rejection(x_list, bad_event_dict)[1]/trials
        iterations.append(avg_iterations)
    return side_lengths, iterations

def exponential_func(x, a, b, c):
    return a*np.exp(-b*x)+c

def regression_analysis(m_low, m_high, n_low, n_high, skip=1, trials = 50, degree = 3, p=0.5):
    """Runs a linear regression to fit a polynomial to the runtimes of GPRS on a square grid.
    Returns a matrix of the runtimes."""
    pass

def exp_fit_analysis(m_low, m_high, n_low, n_high, skip=1, trials=50, p=0.5):
    """Runs an exponential curve fit to the runtimes of GPRS on a square grid.
    b negative implies a logarithmic runtime."""
    for m in range(m_low, m_high+1):
        x, y = grid_results_graph(m, trials, n_low, n_high, skip, p)
        results = curve_fit(exponential_func, x, y)[0]
        print('For m = '+str(m)+", best fit is "+str(results[0])+"e^(-"+str(results[1])+"x)+"+str(results[2]))
        xx = np.linspace(1,30)
        yy = exponential_func(xx, *results)
        plt.plot(x,y,'o', xx, yy)
        plt.show

exp_fit_analysis(13, 13, 1, 30)

