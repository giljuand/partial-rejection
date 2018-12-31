import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation as animation
from partial_rejection import *
from square_grid import *


"""
Simulates general partial rejection sampling on a square grid of dimensions (M,N)
with edge probability P.
"""

N = 100
M = 100
P = 0.45
SPEED = 300

BAD = 2
RES = 1
NOT_RES = 0

grid = np.zeros((M,N))

#the next line is necessary in initialization for some unknown reason
grid[0,0]=2

def square_gprs_iter(n=N, m=M, p=P):
    """Performs general partial rejection sampling given a
    list of the variables and a dictionary mapping bad events to a dictionary mapping variables to values.
    Returns the set of bad events and set of variables to resample every round of Algorithm 6."""
    x_list = create_grid_x_list(n, m, p)
    bad_event_dict = create_grid_bad_event_dict(n,m)
    resample = set(range(len(x_list)))
    bad_events = set()
    assignment = {} #maps each variable to its current assigned value
    iterations = 0

    #constructs var_dict from bad_event_dict
    var_dict = {}
    for A in bad_event_dict:
        for var in bad_event_dict[A]:
            if var not in var_dict:
                var_dict[var] = set()
            var_dict[var].add(A)

    while resample != set():
        yield (resample, bad_events)
        for x in resample:
            assignment[x] = take_sample(x_list[x])
        bad_events = find_bad_events(assignment, bad_event_dict)
        resample = find_resample_set(bad_events, assignment, bad_event_dict, var_dict)
        iterations += 1
    yield (resample, bad_events)

def update(data):
    """
    Takes the variables to be resampled and the bad events
    to create a grid in which Bad, Res\Bad, and not Res are
    labeled by 2,1,0 respectively.
    Returns the matrix figure."""
    newGrid = np.zeros((M, N))
    resample, bad_events = data
    for i in range(M):
        for j in range(N):
            val = N*i+j
            #labels as bad event
            if val in bad_events:
                newGrid[i, j] = BAD
            #labels as res \ bad
            elif (i*(2*N+1)+N+j in resample and i*(2*N+1)+N+j+1 in resample
                    and i*(2*N+1)+j in resample and (i+1)*(2*N+1)+j in resample):
                newGrid[i, j] = RES
            #labels as not resample
            else:
                newGrid[i,j] = NOT_RES
    #update data
    mat.set_data(newGrid)
    return [mat]

#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

# set up animation
fig, ax = plt.subplots()
mat = ax.matshow(grid)
ani = animation.FuncAnimation(fig, update, frames = square_gprs_iter, interval=SPEED,
                              save_count=200, repeat = False, blit=True)
#ani.save('square1.htm', writer = writer)
plt.show()
