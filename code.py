from collections.abc import Set
from math import comb, ceil, floor
import numpy as np
from scipy.optimize import linprog
from collections import Counter
import sys
from pathlib import Path
import re
from time import perf_counter

def file_order(path): # function to read the files in the corresponding benchmark folder
                      # following the listed order
    return [
        int(part) if part.isdigit() else part.casefold()
        for part in re.split(r"(\d+)", str(path))
    ]


def read_binpack(filename): # function to read individual instance txt files and return
                            # the corresponding input parameters

    with open(filename, "r") as f:
        tokens = iter(f.read().split())

    num_sizes = int(next(tokens))
    cap = int(next(tokens))
    sizes = [int(next(tokens)) for _ in range(num_sizes)]

    return [cap, num_sizes, sizes]


def group_sizes(sizes): # group together item sizes, make a width/size and demand list : turn it into a cutting stock problem

    counts = Counter(sizes)        # returns a item_width : demand dictionary, e.g. {33: 4, 42: 2, ...}
    widths = sorted(counts.keys(), reverse=True)
    demands = [counts[w] for w in widths]
    m = len(widths)

    return [m, widths, demands]

# KANTOROVICH MODEL

def first_fit_decreasing(W, sizes): # implements the first-fit-decreasing heuristic to find an upper bound K
                                    # for the number of bins used in the Kantorovich model 

    (m, sizes, demands) = group_sizes(sizes)
    used_bins = 0
    rem_cap_list = []
    for i in range(m):
        for _ in range(demands[i]):
            flag = 0
            for j in range(len(rem_cap_list)):
                if rem_cap_list[j] >= sizes[i]:
                    rem_cap_list[j] -= sizes[i]
                    flag = 1
                    break
                
            if not flag:
                used_bins +=1
                assert W >= sizes[i]
                rem_cap_list.append(W - sizes[i]) 

    return used_bins  # also equal to len(rem_cap_list)


def kantorovich(W, m, sizes, time_limit):

    model_start = perf_counter()

    K = first_fit_decreasing(W, sizes)
    m, sizes, demands = group_sizes(sizes)

    

    # set the objective vector
    c = np.zeros((m * K + K))

    for i in range(-1, -K - 1, -1):
        c[i] = 1

    # set the constraint matrix and vector
    A_ub = np.zeros((m + K, m * K + K))
    b_ub = np.zeros((m + K))

    # fill up the constrain matrix and vector
    for i in range(0, m):
        for j in range(0, K):
            A_ub[i][i * K + j] = -1

        b_ub[i] = - demands[i]

    for i in range(m, m + K):
        for p in range(0, m):
            A_ub[i][p * K + i - m] = sizes[p]

        
        A_ub[i][m * K + i - m] = -W

    # setting up the bounds and integrality for the decision variables
    bounds = [(0, None)] * (m * K) + [(0, 1)] * K
    integrality = np.ones(m * K + K , dtype = int)


    # return the optimum of the ILP using HiGHS solver

    solve_start = perf_counter()

    result = linprog(c, A_ub = A_ub, b_ub = b_ub, bounds = bounds, method = 'highs', integrality = integrality ,
        options={"time_limit": float(time_limit),  # seconds per instance
        "mip_rel_gap": 0.0})

    solve_end = perf_counter()

    result.total_time = solve_end - model_start
    result.ilp_time = solve_end - solve_start

    lp_result = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs", options={"time_limit": float(time_limit)}
)
    if result.success and lp_result.success:
        result.ip_lp_ratio = result.fun / ceil(lp_result.fun)
    else:
        result.ip_lp_ratio = None

    return result



# GILMORE - GOMORY MODEL
def generate_maximal_patterns(W, m, sizes): # generating maximal pattern both preserves the optimal objective and saves on number on varables/columns
                                            # in LP of Gilmore - Gomory

    modified = group_sizes(sizes)
    m = modified[0]
    widths = modified[1]

    # we will count number of patterns by dynamic programming in time O(m * W^2) and space 0(m * W)
    # and then the list of all maximal patterns in O(m^2 * P + m * W * P), where P is number of maximal patterns
    # initiate a (m + 1) X (W + 1) array
    dp = np.zeros((m + 1, W + 1), dtype = object)

    # dp[i][r] represents number of ways to coplete a MAXIMAL pattern using types i,...,m - 1, with r bin-capacity remaining

    # now we fill in the base cases i.e. dp[m][r] which represents all bins are processed with still remaining capacity r
    for r in range(0, W + 1):

        if r < widths[-1]: # if the remaining capacity is less than the smallest width, it is indeed a maximal pattern
                           # and we count 1 for it
            dp[m][r] = 1   
        else:              # if the remaining capacity is equal to or more than the smallest width, it is NOT a maximal pattern
                           # and we count 0 for it
            dp[m][r] = 0


    # table recursion
    for i in range (m - 1, -1, -1):
        for r in range(0, W + 1):
            num = r // widths[i]
            dp[i][r] = sum([dp[i + 1][r - j * widths[i]] for j in range(0, num + 1)])

    # now we compute the list of patterns using dp(i.e. the counts)
    # the tactic is this: process items from 0 to m-1. If say [p_1, ..., p_z] is the (uncompleted) list of patterns at the (i - 1)th stage
    # for the ith stage see which q(# number of items) are feasible for each p_k i.e. doesn't overshoot remaining capacity : q * widths[i] <= remaining_capacity, and,
    # corresponds to  (extension to) a maximal pattern via dp[i + 1][remaining_capacity] > 0 ), then append every such q to p_k to extend p_k

    patterns = [[q] for q in range(0, (W // widths[0]) + 1) if dp[1][W - q * widths[0]] > 0] # I fill out for i = 0 first
    for i in range(1, m): # process from item 1 to m - 1
        extension = []

        for j in patterns: # go through the present uncompleted patterns
            sum_j = 0
            for k in range(0, i):
                sum_j += j[k] * widths[k]  # compute their sum

            remaining_capacity = W - sum_j  # compute remaining capacity

            num_item = 0  # check feasible num_item(s)
            while remaining_capacity >= 0: # num_item number of item i doesn't overshoot remaining capacity

                if dp[i + 1][remaining_capacity] > 0: # num_item number of item i extends to a maximal pattern
                    extension.append(j + [num_item]) 

                num_item += 1
                remaining_capacity -= widths[i] # decrease remaining capacity

        patterns = extension # set the new patterns to be the extension of the old one
            


    return [dp[0][W], patterns] # return number of maximal patterns and all the maximal patterns. NOTE: dp[0][W] = len(patterns)

def ggomory(W, m, sizes, time_limit):

    model_start = perf_counter()

    m, sizes, demands = group_sizes(sizes)


    [num_patterns, patterns] = generate_maximal_patterns(W, m, sizes)
    

    c = np.ones((num_patterns))

    A_ub = np.zeros((m, num_patterns))
    b_ub = np.zeros((m))

    for i in range(0, m):
        for p in range(0, num_patterns):
            A_ub[i][p] = - patterns[p][i]
        b_ub[i] = -demands[i]


    integrality = np.ones(num_patterns , dtype = int)

    solve_start = perf_counter()

    result = linprog(c, A_ub = A_ub, b_ub = b_ub, bounds = (0, None), method = 'highs', integrality = integrality ,
            options={"time_limit": float(time_limit),  # seconds per instance
            "mip_rel_gap": 0.0})

    solve_end = perf_counter()
    
    result.total_time = solve_end - model_start
    result.ilp_time = solve_end - solve_start
    
    lp_result = linprog(c, A_ub = A_ub, b_ub = b_ub, bounds = (0,None), method = "highs", options={"time_limit": float(time_limit)})

    if result.success and lp_result.success:
        result.ip_lp_ratio = result.fun / ceil(lp_result.fun)
    else:
        result.ip_lp_ratio = None
    
    return result
    

# ARCFLOW MODEL
def arc_reduction1(W, out_adj, in_adj):
    
    for i in range(1, W):
        temp1 = out_adj[i].copy()
        for j in temp1:
            if not any(i - k >= j - i for k in in_adj[i]):
                out_adj[i].remove(j)
                in_adj[j].remove(i)
    return (out_adj, in_adj)


def arc_reduction2(W, m, out_adj, in_adj, sizes, demands):

    for i in range(m):
        if sizes[i] == 1:
            continue
        anchor_points = [0]
        for j in range(1, W) :
            if any(j - k > sizes[i] for k in in_adj[j]):
                anchor_points.append(j)
        allowed_starting_points = set()
        for p in anchor_points:
            for z in range(demands[i]):
                if p + z * sizes[i] <= W:
                    allowed_starting_points.add(p + z * sizes[i])

        for node in range(W + 1):
            if node + sizes[i] <= W:
                target = node + sizes[i]
                if (node not in allowed_starting_points) and (target in out_adj[node]):
                
                    out_adj[node].remove(target)
                    in_adj[target].remove(node)




    
    return (out_adj, in_adj)

def arcflow(W, m, sizes, time_limit):

    model_start = perf_counter()

    m, sizes, demands= group_sizes(sizes)

    out_adj = {u: set() for u in range(W + 1)}
    in_adj = {u: set() for u in range(W + 1)}


    # add the arcs
    for i in range(W + 1): # go over vertices 0,..,W and at vertex i add all possible in/out-neighbors corresponding to
                           # different sizes
        for j in range(m):
            if (i + sizes[j] <= W): 
                out_adj[i].add(i + sizes[j])
                in_adj[i + sizes[j]].add(i)

    if sizes[-1] <= W :
        for i in range(sizes[-1], W):
            out_adj[i].add(i + 1)
            in_adj[i + 1]. add(i)

    (out_adj, in_adj) = arc_reduction1(W, out_adj, in_adj) # apply reduction rule 1
    (out_adj, in_adj) = arc_reduction2(W, m, out_adj, in_adj, sizes, demands) # apply reduction rule 2

    out_adj_list = [ list(out_adj[i]) for i in range(W + 1)]
    cumulative_degree_sum_list = [0] * (W + 1)
    for i in range(1, W + 1):
        cumulative_degree_sum_list[i] = cumulative_degree_sum_list[i - 1] + len(out_adj_list[i - 1])

    # set up the objective vector
    c = np.zeros((cumulative_degree_sum_list[W] + len(out_adj_list[W])) + 1)
    c[-1] = 1

    # set up the equality constraint matrix and vector
    A_eq = np.zeros((W + 1, cumulative_degree_sum_list[W] + len(out_adj_list[W]) + 1))
    b_eq = np.zeros((W + 1))

    # fill up the equality constraint matrix and vector
    for i in range(0, W + 1):
        if i == 0:
            for j in out_adj[0]:
                A_eq[0][cumulative_degree_sum_list[0] + out_adj_list[0].index(j)] = 1
            A_eq[0][-1] = -1

        elif i == W :
            for j in in_adj[W]:
                A_eq[W][cumulative_degree_sum_list[j] + out_adj_list[j].index(W)] = -1
            A_eq[W][-1] = 1

        else :
            for j in out_adj[i]:
                A_eq[i][cumulative_degree_sum_list[i] + out_adj_list[i].index(j)] = -1
            for j in in_adj[i]:
                A_eq[i][cumulative_degree_sum_list[j] + out_adj_list[j].index(i)] = 1

    # set up the inequality constraint matrix and vector
    A_ub = np.zeros((m, cumulative_degree_sum_list[W] + len(out_adj_list[W]) + 1))
    b_ub = np.zeros((m))

    # fill up the inequality constraint matrix and vector
    for i in range(0, m):
        for j in range(0, W + 1):
            if (j + sizes[i] in out_adj[j]):
                A_ub[i][cumulative_degree_sum_list[j] + out_adj_list[j].index(j + sizes[i])] = -1
        b_ub[i] = -demands[i]

    integrality = np.ones((cumulative_degree_sum_list[W] + len(out_adj_list[W]) + 1), dtype = int)

    solve_start = perf_counter()
    result = linprog(c, A_ub = A_ub, b_ub = b_ub, A_eq = A_eq, b_eq = b_eq, bounds = (0, None), method = 'highs', integrality = integrality,
                 options = {"time_limit": float(time_limit),  # seconds per instance
                "mip_rel_gap": 0.0})
    solve_end = perf_counter()
    
    result.total_time = solve_end - model_start
    result.ilp_time = solve_end - solve_start
    
    lp_result = linprog(c, A_ub = A_ub, b_ub = b_ub, A_eq = A_eq, b_eq = b_eq, bounds = (0,None), method = "highs", options = {"time_limit": float(time_limit)})
        
    if result.success and lp_result.success:
        result.ip_lp_ratio = result.fun / ceil(lp_result.fun)
    else:
        result.ip_lp_ratio = None
    
    return result


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("USAGE: python code.py <folder_name> <model_name> <time_limit_per_instance_in_seconds>")

    folder = sys.argv[1]
    model = sys.argv[2]
    time_limit = sys.argv[3]
    
    files = sorted(Path(folder).rglob("*.txt"), key = file_order)
    num_instances = len(files)
    instances_solved = 0

    avg_total_time = 0
    avg_ilp_time = 0
    avg_bound_strenth = 0

    
    if model == 'kantorovich':
        for i in range(num_instances):
            W, m, sizes = read_binpack(files[i])

            result = kantorovich(W, m, sizes, time_limit)
            print(f"Instance {i + 1}/{num_instances}: "f"{result.message}", flush=True)
            print("Construction + ILP solving:", result.total_time, "seconds")
            avg_total_time += result.total_time
            print("Only ILP solving:", result.ilp_time, "seconds")
            avg_ilp_time += result.ilp_time
            

            if result.success:
                instances_solved += 1
                print(f"ILP optimum: {round(result.fun)}")
                print("IP/LP:", result.ip_lp_ratio)
                avg_bound_strenth += result.ip_lp_ratio

    elif model == 'gilmore-gomory':
        for i in range(num_instances):
            W, m, sizes = read_binpack(files[i])
        
            result = ggomory(W, m, sizes, time_limit)
            print(f"Instance {i + 1}/{num_instances}: "f"{result.message}", flush=True)
            print("Construction + ILP solving:", result.total_time, "seconds")
            avg_total_time += result.total_time
            print("Only ILP solving:", result.ilp_time, "seconds")
            avg_ilp_time += result.ilp_time
        
            if result.success:
                instances_solved += 1
                print(f"ILP optimum: {round(result.fun)}")
                print("IP/LP:", result.ip_lp_ratio)
                avg_bound_strenth += result.ip_lp_ratio

    elif model == 'arcflow':
        for i in range(0, num_instances):
            W, m, sizes = read_binpack(files[i])
            result = arcflow(W, m, sizes, time_limit)
            print(f"Instance {i + 1}/{num_instances}: "f"{result.message}", flush=True)
            print("Construction + ILP solving:", result.total_time, "seconds")
            avg_total_time += result.total_time
            print("Only ILP solving:", result.ilp_time, "seconds")
            avg_ilp_time += result.ilp_time
                    
                    
            if result.success:
                instances_solved += 1
                print(f"ILP optimum: {round(result.fun)}")
                print("IP/LP:", result.ip_lp_ratio)
                avg_bound_strenth += result.ip_lp_ratio

    else:
        sys.exit("USAGE: only allowed models are kantorovich, gilmore-gomory, arcflow")

    print(f"Total instances solved to optimality: {instances_solved} out of {num_instances}")
    print(f"Average total time: {avg_total_time / num_instances} seconds")
    print(f"Average ILP time: {avg_ilp_time / num_instances} seconds")
    if instances_solved == 0:
        print(f"Average bound strength: {None}")
    else:
        print(f"Average bound strength: {avg_bound_strenth / instances_solved}")

        


