import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation as animation
from partial_rejection import *
from honeycomb import *
from matplotlib import collections as mc


"""
Simulates general partial rejection sampling on a square grid of dimensions (M,N)
with edge probability P.

Includes a graph demonstrating size of resample set and bad events over time.
"""

N = 50
M = 50
P = .46
SPEED = 300


def create_edges(n, m, assignment):
    if assignment == {}:
        return []

    lines = []
    for i in range((3*n+1)*m+n):
        if assignment[i]==True:
            r = i//(3*n+1) #row
            p = i%(3*n+1) #position
            if p<n:
                lines.append([((m-r)/2+p, (m-r)*3**0.5/2), ((m-r)/2+p+1, (m-r)*3**0.5/2)])
            else: #in vertical range
                lines.append([((m-r)/2+(p-n)//2, (m-r)*3**0.5/2), ((m-r-1)/2+(p-n+1)//2, (m-r-1)*3**0.5/2)])
    return lines


def honeycomb_gprs_iter(n=N, m=M, p=P):
    """Performs general partial rejection sampling given a
    list of the variables and a dictionary mapping bad events to a dictionary mapping variables to values.
    Returns the set of bad events and set of variables to resample every round of Algorithm 6."""
    x_list = create_honeycomb_x_list(n, m, p)
    bad_event_dict = create_honeycomb_bad_event_dict(n,m)
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
edges.set_xlim(-1, M/2+N+1)
edges.set_ylim(-1, 3**0.5/2*M+1)
plt.gca().set_aspect('equal', adjustable='box')


ani = animation.FuncAnimation(fig, update, frames = honeycomb_gprs_iter, interval=SPEED,
        save_count=5, repeat = False, blit=False) #note: blit must be set to False
plt.show()
