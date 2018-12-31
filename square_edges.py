import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation as animation
from partial_rejection import *
from square_grid import *
from matplotlib import collections as mc


"""
Simulates general partial rejection sampling on a square grid of dimensions (M,N)
with edge probability P.

Includes a graph demonstrating size of resample set and bad events over time.
"""

N = 50
M = 50
P = 0.1
SPEED = 300


def create_edges(n, m, assignment):
    if assignment == {}:
        return []

    lines = []
    horizontal = list(range(n))
    for i in range((2*n+1)*m+n):
        if assignment[i]==True:
            r = i//(2*n+1) #row
            p = i%(2*n+1) #position
            if p in horizontal:
                lines.append([(p, m-r), (p+1, m-r)])
            else: #in vertical range
                lines.append([(p-n, m-r), (p-n, m-r-1)])
    return lines


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
        yield (resample, bad_events, iterations, assignment)
        for x in resample:
            assignment[x] = take_sample(x_list[x])
        bad_events = find_bad_events(assignment, bad_event_dict)
        resample = find_resample_set(bad_events, assignment, bad_event_dict, var_dict)
        iterations += 1
    yield (resample, bad_events, iterations, assignment)

def update(data):
    """
    Takes the variables to be resampled and the bad events
    to create a grid in which Bad, Res\Bad, and not Res are
    labeled by 2,1,0 respectively.
    Returns the matrix figure."""
    resample, bad_events, iterations, assignment = data
    new_edges = create_edges(N, M, assignment)
    new_lc = mc.LineCollection(new_edges)
    
    #update data
    lc.set_segments(new_edges)

# set up animation
fig, axes = plt.subplots()

lc = mc.LineCollection([])
edges = axes
edges.add_collection(lc)
edges.set_xlim(-1, N+1)
edges.set_ylim(-1, M+1)
plt.gca().set_aspect('equal', adjustable='box')


ani = animation.FuncAnimation(fig, update, frames = square_gprs_iter, interval=SPEED,
        save_count=5, repeat = False, blit=False) #note: blit must be set to False
plt.show()
