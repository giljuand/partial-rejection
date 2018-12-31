from partial_rejection import *

# CREATION OF GRAPHS (in x_list form)
# Note:  roots are excluded since no associated variable

def create_cycle(n):
    result=[]
    for i in range(n-1):
        result.append([((i-1)%n, 1),((i+1)%n, 1)])
    return result

def create_weighted_cycle(n, m=5):
    """Creates a cycle with one edge (0,1) that is m times heavier than the others.
    Note:  n must be 3 or larger."""
    result = []
    result.append([(n-1, 1), (1, m)])
    result.append([(0, m), (2, 1)])
    for i in range(2,n-1):
        result.append([((i-1)%n, 1),((i+1)%n, 1)])
    return result

def create_complete_graph(n):
    result=[]
    for i in range(n-1):
        result.append([(j, 1) for j in range(n) if j!=i])
    return result

def create_weighted_complete_graph(n, m=100):
    result = [[] for i in range(n-1)] # n-1 variables (excludes root)
    for i in range(n-1): #goes through variables
        for j in range(i): #goes through variables of smaller index
            weight = random.randint(1,m)
            result[i].append((j, weight))
            result[j].append((i, weight))
        root_weight = random.randint(1,m)
        result[i].append((n-1, root_weight))
    return result

def create_sparse_graph(n, d=0.5, m=10):
    """Creates a graph with a guaranteed cycle of weight 1 edges around the outside, where
    every other edge has probability d of appearing.
    Cycles have random weights between 1 and m."""
    result = create_cycle(n)
    for i in range(n-1): #goes through variables
        for j in range(i): #goes through variables of smaller index
            if i != (j-1)%n and i != (j+1)%n and random.random()<d: #d chance of creating edge
                weight = random.randint(1,m)
                result[i].append((j, weight))
                result[j].append((i, weight))
        if i != n-2 and i != 0 and random.random()<d: #connects to root if not yet connected
            root_weight = random.randint(1,m)
            result[i].append((n-1, root_weight))
    return result


# FINDING CYCLES

def find_cycles(assignment):
    """Given an assignment of vertices to a neighboring vertex, returns a list of variables
    in a cycle and the number of cycles located."""
    vertices_left = list(assignment.keys()) #all nodes except for root
    cycle_vertices = []
    num_cycles = 0
    while vertices_left != []:
        v = vertices_left[0]
        path = []
        cycle_found = False
        while not cycle_found:
            if v in path: #cycle detected
                cycle_vertices += path[path.index(v):] #adds all vertices after v inclusive
                print(path[path.index(v):])
                cycle_found = True
                num_cycles += 1
            elif v not in vertices_left: #no cycle but end of path, e.g. hits root
                break
            else: #no cycle yet
                vertices_left.remove(v)
                path.append(v)
                v = assignment[v]
    return (cycle_vertices, num_cycles)


# MAXIMUM SPANNING TREES

def create_graph_array(x_list, n):
    """Creates an array representing the size n graph given by x_list.  Used to calculate minimal spanning tree."""
    result = [[0 for i in range(n)] for j in range(n)]
    for i in range(len(x_list)):
        for j in x_list[i]:
            result[i][j[0]]=j[1] #j[1] is the weight of the edge between i and j[0]
            result[j[0]][i]=j[1] #keeps results symmetric, since root can be problematic
    return result

def create_unweighted_graph(x_list):
    """Returns the unweighted version of x_list."""
    result = []
    for var in x_list: #goes through each variable
        uw_var = []
        for value in var: #goes through each value the variable can take
            if value[1]!=0: #assigns them the same weight, unless weight is zero
                uw_var.append((value[0], 1))
            else:
                uw_var.append((value[0], 0))
        result.append(uw_var)
    return result

def find_weight(assignment, g=None):
    """Calculates the weight of a spanning tree on g with certain edge assignment, either dict or csr_matrix.
    g is the array representation created by create_graph_array, necessary if assignment is dict."""
    total = 0
    if type(assignment)==dict: #assignment from PRS
        for var in assignment.keys(): 
            total+=g[var][assignment[var]]
    elif type(assignment)==csr_matrix: #assignment from minimum_spanning_tree
        arr = assignment.toarray().astype(int) #turns assignment into array
        total = sum([sum(row) for row in arr])
    return total

def maximum_spanning_tree(g):
    """Returns the weight of the maximum spanning tree of graph g."""
    g1 = [[-val for val in row] for row in g] #negate edge weights
    return - find_weight(minimum_spanning_tree(g1))



#PLOTTING RESULTS
def spanning_tree_runtime_plot(create_graph, n):
    """Creates a plot analyzing the runtime of spanning tree PRS as size of graph increases."""
    TRIALS = 300
    x=[]
    y=[]
    for i in range(1, n):
        x.append(i)
        y.append(expected_bad_events(create_graph(i), find_cycles, TRIALS)) #expected number for graph size i
    plt.plot(x,y, 'ro') #red dots
    plt.xlabel("Size of graph")
    plt.ylabel("Average # bad events over "+str(TRIALS)+" trials")
    plt.title("PRS for complete graph spanning tree")
    plt.show()

def PRS_max_spanning_tree_comparison(create_graph, t, n):
    """Creates a plot comparing the minimum spanning tree of a graph with
    the minimum spanning tree created over t trials of weighted PRS,
    the average spanning tree created over t trials of weighted PRS,
    the average spanning tree created over t trials of unweighted PRS,
    and the minimum spanning tree created over t trials of unweighted PRS,
    as the size of the graph increases up to size n.
    Note:  must be weighted"""

    ratio = []
    for i in range(30, n):
        x_list = create_graph(i)
        uwx_list = create_unweighted_graph(x_list)
        g = create_graph_array(x_list, i)
        uw_PRS_weights = [] #the weights of max spanning trees using unweighted PRS
        w_PRS_weights = [] #the weights of max spanning trees using weighted PRS

        #runs the t trials of PRS
        for j in range(t):
            uw_tree = partial_rejection(uwx_list, find_cycles)[0]
            w_tree = partial_rejection(x_list, find_cycles)[0]
            uw_PRS_weights.append(find_weight(uw_tree, g))
            w_PRS_weights.append(find_weight(w_tree, g))
        avg_uw=sum(uw_PRS_weights)/t
        avg_w = sum(w_PRS_weights)/t
        max_uw = max(uw_PRS_weights)
        max_w = max(w_PRS_weights)
        max_tree = maximum_spanning_tree(g)
        plt.plot(i, avg_uw, 'go', i, avg_w, 'bo', i, max_uw, 'g+', i, max_w, 'b+', i, max_tree, 'ro')
        ratio.append(max_tree/max_w)

    #creation of legend
    red_patch = mpatches.Patch(color = 'red', label = 'Maximum tree')
    green_patch = mpatches.Patch(color = 'green', label = 'unweighted PRS')
    blue_patch = mpatches.Patch(color = 'blue', label = 'weighted PRS')
    plt.legend(handles=[red_patch, green_patch, blue_patch])

    plt.xlabel('Size of randomly-weighted complete graph')
    plt.ylabel('Weight of spanning tree obtained')
    plt.title('Weighted spanning trees in PRS')
    print(ratio)
    print(sum(ratio)/len(ratio))
    plt.show()

def spanning_tree_weight_distribution(create_graph, t, n):
    weights = []
    x_list = create_graph(n)
    g = create_graph_array(x_list, n)
    for i in range(t):
        tree = partial_rejection(x_list, find_cycles)[0]
        weights.append(find_weight(tree, g))
    plt.hist(weights)
    print("max tree:", maximum_spanning_tree(g))
    print('avg tree:', sum(weights)/t)
    plt.show()


def one_minus_epsilon_tree(create_graph, epsilon, n, t=1):
    x = [] #size of graph
    y = [] #iterations until (1-epsilon)MST achieved
    for i in range(5, n):
        print(i)
        iteration_vals = []
        for j in range(t):
            x_list = create_graph(i)
            g = create_graph_array(x_list, i)
            mst = maximum_spanning_tree(g)
            w = -float("inf")
            iterations=0
            while w<(1-epsilon)*mst:
                w = find_weight(partial_rejection(x_list, find_cycles)[0], g)
                iterations+=1
            iteration_vals.append(iterations)
        x.append(i)
        y.append(sum(iteration_vals)/t)
    print(x,y)
    print(len(x))
    plt.plot(x,y,'ro')
    plt.xlabel("Size of complete graph")
    plt.ylabel("Average iterations until weight > " + str(1-epsilon) + " * MST over " + str(t)+" trials")
    plt.show()
