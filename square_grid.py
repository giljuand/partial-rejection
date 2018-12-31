from partial_rejection import *
import time


def create_grid_x_list(n,m, p=0.5):
    """Creates the variables for a square grid with dimensions
    n and m.  w is how much more likely an edge is to not be there."""
    result = []
    for i in range((2*n+1)*m+n):
        result.append([(True,p),(False,1-p)])
    return result

def create_grid_bad_event_dict(n,m):
    result = {}
    for i in range(n*m): #bad events
        r = i//n
        c = i % n
        var_dict = {r*(2*n+1)+n+c: True, 
                r*(2*n+1)+n+c+1: True,
                r*(2*n+1)+c:  True,
                (r+1)*(2*n+1)+c: True}
        result[i]=var_dict
    return result

def grid_results_graph(trials, low, high, skip=1, p=0.5) :
    side_lengths = []
    iterations = []
    done = False #cuts off when slows down
    for i in range(low, high, skip):
        print('Size of graph: ', i)
        avg_iterations = 0
        start_time = time.time()
        for j in range(trials):
            #print('Trial: ', j)
            x_list = create_grid_x_list(i, i, p)
            bad_event_dict = create_grid_bad_event_dict(i, i)
            avg_iterations += general_partial_rejection(x_list, bad_event_dict)[1]/trials
            if time.time()-start_time>150: #time until cut off
                done = True
                break
        if done: #simulations taking too long, return graph
            break
        iterations.append(avg_iterations)
        side_lengths.append(i)
    plt.xlabel('Square grid of size n x n, p='+str(p))
    plt.ylabel('Average rounds of G.P.R.S over '+str(trials)+' trials')
    plt.plot(side_lengths, iterations, 'ro')
    plt.savefig('PRS {0}.png'.format(p))
    plt.clf()

"""
i = 0.25
tested_probs = []
x = True
while x:
    j = 0.25
    while j<=0.5:
        if j not in tested_probs:
            print(j)
            grid_results_graph(200, 5, 80, 1, j)
            tested_probs.append(j)
        j+=i
    x = False
    i = i/2
"""
