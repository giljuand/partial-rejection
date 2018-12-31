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

N = 10
M = 10
P = 0.47
SPEED = 1000

BAD = 2
RES = 1
NOT_RES = 0

#initialize grid
grid = np.zeros((M,N))

#initialize values for graph
xs=[]
num_res=[]
num_bad=[]

#the next line is necessary in initialization for some unknown reason
grid[0,0]=2


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
    newGrid = np.zeros((M, N))
    resample, bad_events, iterations, assignment = data
    new_edges = create_edges(N, M, assignment)
    new_lc = mc.LineCollection(new_edges)
    bad_size = 0
    res_size = 0
    for i in range(M):
        for j in range(N):
            val = N*i+j
            #labels as bad event
            if val in bad_events:
                newGrid[i, j] = BAD
                bad_size+=1
            #labels as res \ bad
            elif (i*(2*N+1)+N+j in resample and i*(2*N+1)+N+j+1 in resample
                    and i*(2*N+1)+j in resample and (i+1)*(2*N+1)+j in resample):
                newGrid[i, j] = RES
                res_size+=1
            #labels as not resample
            else:
                newGrid[i,j] = NOT_RES
    #update data
    if iterations!=0: #updates lists of data
        xs.append(iterations)
        num_bad.append(bad_size)
        num_res.append(res_size+bad_size)
    r.set_data(xs, num_res) #loads data for resample set
    b.set_data(xs, num_bad) #loads data for bad events
    sg.relim() #clears limiters in lower subplot
    sg.autoscale_view(True) #scales the axes

    lc.set_segments(new_edges)

# set up animation
fig, axes = plt.subplots(2,1)

lc = mc.LineCollection([])
edges = axes[0]
edges.add_collection(lc)
edges.set_xlim(-1, N+1)
edges.set_ylim(-1, M+1)
#edges.axis('equal')

sg = axes[1] #size graph
r, = sg.plot([],[], 'b-') #initializes line object for resample set
b, = sg.plot([],[], 'r-') #initializes line object for bad events
sg.set_title('Size of sets over iterations')
sg.set_xlabel('Iterations')
sg.set_ylabel('Size of set')
sg.legend([r,b],['Resample set', 'Bad events'])

ani = animation.FuncAnimation(fig, update, frames = square_gprs_iter, interval=SPEED,
        save_count=5, repeat = False, blit=False) #note: blit must be set to False
plt.show()
