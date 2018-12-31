from partial_rejection import *
import math
import time

def create_honeycomb_x_list(n,m, p=0.5):
    """Creates the variables for a honeycomb graph with dimensions
    n and m.  w is how much more likely an edge is to not be there."""
    result = []
    for i in range((3*n+1)*m+n):
        result.append([(True,p),(False,1-p)]) #equal chance 
    return result

def create_honeycomb_bad_event_dict(n,m):
    result = {}
    for i in range(2*n*m):
        r = i//(2*n)
        c = i % (2*n)
        var_dict = {r*(3*n+1)+n+c: True, r*(3*n+1)+n+c+1: True} #side edges
        if c%2==1: #top edge
            var_dict[r*(3*n+1)+c//2]=True
        else: #bottom edge
            var_dict[(r+1)*(3*n+1)+c//2]=True
        result[i]=var_dict
    return result

def honeycomb_results_graph(trials, low, high, skip=1, p=0.5) :
    side_lengths = []
    iterations = []
    done = False
    for i in range(low, high, skip):
        print('Size of graph: ', i)
        avg_iterations = 0
        start_time = time.time()
        for j in range(trials):
            #print('Trial: ', j)
            x_list = create_honeycomb_x_list(i, i, p)
            bad_event_dict = create_honeycomb_bad_event_dict(i, i)
            avg_iterations += general_partial_rejection(x_list, bad_event_dict)[1]/trials
            if time.time()-start_time>120: #time until cut off
                done = True
                break
        if done: #simulations taking too long, return graph
            break
        iterations.append(avg_iterations)
        side_lengths.append(i)
    plt.xlabel('Triangle grid of size n x n, p='+str(p))
    plt.ylabel('Average rounds of G.P.R.S over'+str(trials)+' trials')
    plt.plot(side_lengths, iterations, 'ro')
    plt.savefig('PRS {0}.png'.format(p))
    plt.clf() #clears figure


#honeycomb_results_graph(200,5,80,1,0.25)

"""
i = 0.001
tested_probs = [0.47, 0.470]
while True:
    j = 0.468
    while j<=0.472:
        if j not in tested_probs:
            print(j)
            honeycomb_results_graph(200, 5, 150, 1, j)
            tested_probs.append(j)
        j+=i
    i = i/2
"""
