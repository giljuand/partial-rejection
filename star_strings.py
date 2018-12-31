from partial_rejection import *

def create_star_string_x_list(n, alphabet):
    result = []
    for x in range(n):
        result.append([(i,1) for i in range(alphabet)])
    return result

def create_star_string_bad_event_dict(n, w):
    """w is a string (with stars) represented as a list."""
    result = {}
    for i in range(n-len(w)+1):
        var_values = {}
        for j in range(len(w)):
            if w[j] != '*':
                var_values[i+j] = int(w[j])
        result[i]=var_values
    return result

def w_free_strings_graph(w, alphabet, trials, low, high, skip=1) :
    """Runs GPRS to sample w_free strings and returns a graph showing expected rounds over the number of trials.
    w is the substring to be avoided, in the form of a list or a string if alphabet <= 10.
    alphabet is the number of characters in the alphabet (e.g. 10 will use digits 0-9)
    low and high are the string lengths to be sampled."""
    lengths = []
    iterations = []
    for i in range(low, high, skip):
        avg_iterations = 0
        lengths.append(i)
        for j in range(trials):
            x_list = create_star_string_x_list(i, alphabet)
            bad_event_dict = create_star_string_bad_event_dict(i, w)
            avg_iterations += general_partial_rejection(x_list, bad_event_dict)[1]/trials
        iterations.append(avg_iterations)
    plt.plot(lengths, iterations, 'ro')
    plt.show()

w_free_strings_graph('11111', 3, 2000, 30, 100)
