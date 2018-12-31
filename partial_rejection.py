import random
import numpy
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
import copy
import matplotlib.patches as mpatches


"""
This code implements partial rejection sampling and general partial rejection sampling.

To apply GPRS to a certain problem, you must define the x_list and the bad_event_dict.

x_list is a list of the variables and the values that the variables can take.

The ith index of x_list corresponds to the variable x_i.

The variables are represented by lists of the values that it can take.

The values are represented by tuples where the zeroeth index is the value and
the first index is the probability or "weight" given to that value.

For example, a variable x_i that takes heads with probability 0.5 and
tails with probability 0.5 might be represented as [('H', 0.5), ('T', 0.5)]
or [('H', 1), ('T', 1)] or [('H', 100), ('T', 100)]. These representations are equivalent.

A biased coin x_i is represented as [('H', 4), ('T', 1)].  Here, we get tails with probability 0.2.


bad_event_dict is a dictionary where the keys are the bad events.  The key i corresponds
to bad event A_i.

The key values are also dictionaries.  In these dictionaries, the keys are the variables
(where j represents variables x_j) and the key values are the values that must be taken
by the variables for the bad event to occur.


In a problem with one coin flip in which we wish to avoid heads, we would have the following:
    x_list = [[('H', 1), ('T', 1)]]
    bad_event_dict = {0: {0: 'H'}}

"""






# PARTIAL REJECTION SAMPLING

def take_sample(x):
    """x is a variable of form e.g.[('H', 1),('T', 1)].
    We pick result with probability directly proportional to weight.
    Outputs a random result."""
    z = 0 #sum of weights
    for i in range(len(x)):
        z+=x[i][1]
    random_num = random.uniform(0,z)
    for i in range(len(x)): #determines what value is associated with the random number
        if random_num <= x[i][1]:
            return x[i][0]
        random_num -= x[i][1]

def partial_rejection(x_list, find_bad_var):
    """Conducts partial rejection sampling.

    x_list is a list in which each index represents a variable with probability distribution based
    on list at that index.  e.g. [[('H', 0.5),('T',0.5)]] is a list with a single fair coin variable.

    find_bad_var takes a sample and should return the variables to be resampled and number of bad events.
    
    Returns a tuple containing the assignment, number of iterations, and number of bad events popped."""

    resample = list(range(len(x_list)))
    assignment = {} #maps each variable to its current assigned value
    iterations = 0
    bad_events = 0
    while resample !=[]:
        print("new round!")
        for x in resample:
            assignment[x] = take_sample(x_list[x])
        bad_var_info = find_bad_var(assignment)
        resample = bad_var_info[0]
        bad_events += bad_var_info[1]
        iterations += 1
    return (assignment, iterations, bad_events)


def find_border(R, N, partial_assignment, variables_to_search, bad_event_dict, var_dict):
    """Given a set of bad events R, a set of bad events N not to be resampled,
    a set of the variables to search for neighbors, and the two remaining dictionaries,
    returns a set of the border bad events to also be resampled.
    Modifies N."""
    unlabeled_neighbor = set() #keeps track of possible border events discovered
    for var in variables_to_search:
        for A in var_dict[var]: #goes through bad events dependent on variable
            if A not in R and A not in N:  #unlabeled neighbor
                if bad_event_dict[A][var]!=partial_assignment[var]: #if disjoint w/ partial assignment
                    N.add(A) #not resample
                    unlabeled_neighbor.discard(A) #A no longer possible border event
                else: #still compatible with partial assignment
                    unlabeled_neighbor.add(A) #A might be in dR
    return unlabeled_neighbor

def find_new_variables(partial_assignment, assignment, border, bad_event_dict):
    """Returns the set of variables associated with bad events at the border
    which are not yet in the partial assignment.
    Modifies partial assignment."""
    result = set()
    for A in border:
        for var in bad_event_dict[A]:
            if var not in partial_assignment:
                result.add(var)
                partial_assignment[var]=assignment[var]
    return result

def find_resample_set(R, assignment, bad_event_dict, var_dict):
    """Finds the resample set used for general partial rejection sampling.
    R is a set of the bad events identified in an assignment.
    assignment is a dictionary mapping variables to assigned values.
    bad_event_dict maps the bad events to a dictionary of variables mapped to their values.
    var_dict maps the variables to a set of the bad events that depend on that variable.
    Returns a set of variables to be resampled."""

    resample = R.copy()
    not_resample = set()

    #constructs partial assignment restricted to R
    partial_assignment = {} #to be returned at end
    for A in R:
        for var in bad_event_dict[A]:
            partial_assignment[var] = assignment[var]

    border = find_border(resample, not_resample, partial_assignment,
                            list(partial_assignment.keys()), bad_event_dict, var_dict)
    resample.update(border) #adds border events to resample set

    while border!=set():
        new_variables = find_new_variables(partial_assignment, assignment, border, bad_event_dict)
        border = find_border(resample, not_resample, partial_assignment, new_variables, bad_event_dict, var_dict)
        resample.update(border)
    return set(partial_assignment.keys())

def find_bad_events(assignment, bad_event_dict):
    """Given an assignment and the bad_event_dict, returns a set of the 
    bad events that occur over that assignment."""
    result = set()
    for A in bad_event_dict:
        if all(bad_event_dict[A][var]==assignment[var] for var in bad_event_dict[A]):
            result.add(A)
    return result

def general_partial_rejection(x_list, bad_event_dict):
    """Performs general partial rejection sampling given a
    list of the variables and a dictionary mapping bad events to a dictionary mapping variables to values."""
    resample = set(range(len(x_list)))
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
        #print(iterations)
        #print('Resample: ', resample)
        #print('Assignment: ', assignment)
        for x in resample:
            assignment[x] = take_sample(x_list[x])
        bad_events = find_bad_events(assignment, bad_event_dict)
        resample = find_resample_set(bad_events, assignment, bad_event_dict, var_dict)
        iterations += 1
    return (assignment, iterations)




# PLOTTING RESULTS

def expected_bad_events(x_list, find_bad_var, trials):
    """Calculates the expected number of bad events resampled, calculated using specified number of trials."""
    total = 0
    for i in range(trials):
        total+=partial_rejection(x_list, find_bad_var)[2]
    return total/trials

def bad_events_histogram(x_list, find_bad_var, trials):
    """Displays a histogram of the number of bad events resampled in specified number of trials."""
    result = []
    for i in range(trials):
        result.append(partial_rejection(x_list, find_bad_var)[2])
    plt.hist(result)
    print(sum(result)/trials)
    plt.show()



