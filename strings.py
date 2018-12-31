from partial_rejection import *

# CREATION OF STRINGS
def create_string_x_list(n, alphabet):
    result = []
    for x in range(n):
        result.append([(i,1) for i in range(alphabet)])
    return result

def create_string_bad_event_dict(n, w):
    result = {}
    for i in range(n-len(w)+1):
        var_values = {}
        for j in range(len(w)):
            var_values[i+j]=int(w[j])
        result[i]=var_values
    return result



# PLOTTING RESULTS
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
            x_list = create_string_x_list(i, alphabet)
            bad_event_dict = create_string_bad_event_dict(i, w)
            avg_iterations += general_partial_rejection(x_list, bad_event_dict)[1]/trials
        iterations.append(avg_iterations)
    plt.plot(lengths, iterations, 'ro')
    plt.show()


def uniformity_test_1(n,alphabet,w, trials):
    """basic test for uniformity of the sampled strings."""
    x_list = create_string_x_list(n, alphabet)
    bad_event_dict = create_string_bad_event_dict(n, w)
    avg=0
    for i in range(trials):
        results=general_partial_rejection(x_list, bad_event_dict)
        sampled_num = ''
        for j in range(n):
            sampled_num += str(results[0][j])
        #print(sampled_num)
        avg += int(sampled_num)/trials
    print(avg)

w_free_strings_graph('11', 2, 20, 8, 8*60, 8)
