# Welcome  to my BPP ILP implementation
## Acknowledgement
This `README.md` file and other stuff has been inspired from [dr. Maxence Delorme](https://www.tilburguniversity.edu/staff/m-delorme)'s repository and his paper:
- M. Delorme, M. Iori, and S. Martello. [BPPLIB: A library for bin packing and cutting stock problems](https://link.springer.com/article/10.1007/s11590-017-1192-z). Optimization Letters, 12(2):235-250, 2018

Other references have been put in the experimental evaluation report available at the end and in-parts(one to be precise) in another section.
## Contents
This repository of the following items:
- The codefile named `code.py` containing the implementation of the models, file-reading/parsing and helper functions among other things as will be explained later
- $4$ benchmark folders named `Falkenauer_U`,`Falkenauer_T`,`many_paterns_unit_demand` and `near_half_widths`
- This `README.md` markdown file
- A `LICENSE` file
## Problem(s) description
**Bin Packing Problem(BPP):** Given $n$ items, each having a integer size $w_j$($j = 1, \ldots, m$) and an unlimited supply of identical bins each with capacity $W$, determine the minimum number of bins required to fit all the items, subject to not surpassing the bins' capacities.  

**Cutting Stock Problem(CSP):** Given $m$ item types, each having an integer weight $w_j$ and an integer demand $d_j$($j=1,\ldots,m$), and an unlimited number of identical bins of integer capacity $W$, determine the minimum number of bins to produce $d_j$ copies of each item type $j$, subject to not surpassing the bins' capacities.
## ILP models for BPP/CSP


### Kantorovich model

Let $W$ be the bin capacity, $m$ the number of distinct item sizes,
$w_i$ the size of item type $i$, and $b_i$ its demand.
Let $K$ be a known upper bound on the number of bins required.

The variable $x_{ik}$ counts the items of type $i$ assigned to bin $k$,
and $y_k$ indicates whether bin $k$ is used.

$$
\begin{aligned}
\min \quad
    & \sum_{k=1}^{K} y_k \\
\text{subject to} \quad
    & \sum_{k=1}^{K} x_{ik} \ge b_i,
    && i=1,\ldots,m, \\
    & \sum_{i=1}^{m} w_i x_{ik} \le W y_k,
    && k=1,\ldots,K, \\
    & y_k \in \{0,1\},
    && k=1,\ldots,K, \\
    & x_{ik} \in \mathbb{Z}_{\geq 0},
    && i=1,\ldots,m,\quad k=1,\ldots,K.
\end{aligned}
$$
### Gilmore–Gomory model with maximal patterns

Let $W$ be the bin capacity, $m$ the number of distinct item sizes,
$w_i$ the size of item type $i$, and $b_i$ its demand.

Let $P_{\max}$ index all maximal cutting patterns.
For each pattern $p$, the coefficient $a_i^p$ counts the items of
type $i$ in that pattern, with $a_i^p \in \mathbb{Z}_{\ge 0}$.

A feasible pattern is maximal if no additional item of any type fits:

$$
0 \le W-\sum_{i=1}^{m} w_i a_i^p
< \min_{1 \le i \le m} w_i.
$$

The variable $\lambda_p$ counts the bins using pattern $p$.

$$
\begin{aligned}
\min \quad
    & \sum_{p \in P_{\max}} \lambda_p \\
\text{subject to} \quad
    & \sum_{p \in P_{\max}} a_i^p \lambda_p \ge b_i,
    && i=1,\ldots,m, \\
    & \lambda_p \in \mathbb{Z}_{\geq 0},
    && p \in P_{\max}.
\end{aligned}
$$

Restricting to maximal patterns preserves the optimal objective:
every feasible pattern can be extended to a maximal one without
using another bin, and the demand constraints allow surplus items.

### Arc-flow model

Use vertices $V=\{0,1,\ldots,W\}$. Add item arcs $(d,d+w_i)$
whenever $d+w_i\le W$, and unit loss arcs $(d,d+1)$ representing
unused capacity.

Let $A'$ contain the forward arcs retained after reduction.
The feedback arc $(W,0)$ is represented separately by $z$.

The variable $x_{de}$ is the flow on arc $(d,e)\in A'$,
while $z$ counts the bins used. Each unit of flow from $0$ to $W$
represents one bin.

$$
\begin{aligned}
\min \quad
    & z \\[4pt]
\text{subject to} \quad
    & \sum_{(0,e)\in A'} x_{0e}-z=0, \\
    & \sum_{(d,v)\in A'} x_{dv}
      -\sum_{(v,e)\in A'} x_{ve}=0,
    && v=1,\ldots,W-1, \\
    & z-\sum_{(d,W)\in A'} x_{dW}=0, \\
    & \sum_{(d,d+w_i)\in A'} x_{d,d+w_i}\ge b_i,
    && i=1,\ldots,m, \\
    & x_{de}\in\mathbb{Z}_{\ge 0},
    && (d,e)\in A', \\
    & z\in\mathbb{Z}_{\ge 0}.
\end{aligned}
$$

**Arc reductions.** Assume $w_1>\cdots>w_m$.

1. **Size ordering:** Keep an item arc $(d,d+w_i)$ only if
   $d=0$ or some retained incoming item arc at $d$ has size
   at least $w_i$.

2. **Remove initial loss:** Delete loss arcs $(d,d+1)$
   with $d<w_m$.

3. **Demand-based reduction:** For type $i$, take anchors $d=0$
   and nodes reached by a strictly larger item arc. Keep type-$i$
   arcs starting at $d+s w_i$, where $s=0,\ldots,b_i-1$,
   provided the entire arc fits within capacity $W$.

In this implementation, reduction 3 skips width-$1$ arcs because
they also serve as loss arcs.
### References:
-  Valerio de Carvalho, J. M., 2002. "[LP models for bin packing and cutting stock problems](https://ideas.repec.org/a/eee/ejores/v141y2002i2p253-273.html)," European Journal of Operational Research, Elsevier, vol. 141(2), pages 253-273, September.
## Benchmarks
All instance files inside the corresponding benchmark folders are given in the following **BPP** format:

- Number of items ($n$) (**NOTE** : $m$ and $n$ are used interchangeably)
- Capacity of the bins ($W$)
- For each item $j$($j = 1,\ldots, m$) 
    - Its size / width($w_i$)
## Usage
Make sure you have Python's **NumPy** and **SciPy** installed. 

To run my code on all the instances on any benchmark folder type the following in the terminal, in the same directory as the code and the folders: 
~~~sh
python code.py <path/to/folder> <model_name> <time_limit_per_instance_in_seconds>
~~~
where `<model_name>` is one of `kantorovich`, `gilmore-gomory` or `arcflow`.

## A walk through my codefile

```python
def file_order(path):
```
This helper function is used to parse the `.txt` files in the benchmark folders in their listed order given the `path` to the folder

~~~python
def read_binpack(filename):
~~~ 
This helper function reads individual instance `.txt` files given the `filename` and returns the corresponding input data parameters to be used by other functions.

~~~python
def first_fit_decreasing(W, sizes):
~~~ 
This helper function implements the first-fit-decreasing heuristic, given `W`(bin-capacity) and `sizes`(item size list) to find an upper bound $K$ for the number of bins used in the Kantorovich model.

~~~python
def kantorovich(W, m, sizes, time_limit):
~~~ 
This functions implements the Kantorovich model given `W`(bin capacity), `m`(number of items), `sizes`(the item size list) and `time_limit`(the time limit in seconds per instance). I use 
`
model_start = perf_counter()
`
at the beginning to track the total runtime including model making. It is matched with `solve_start = perf_counter()` and `solve_end = perf_counter()` to track the ILP solving time. The sizes are grouped via `group_sizes(sizes)` to convert into a **CSP(Cutting Stock Problem)** instance with new `m` and `sizes` and a demands lsit `demands`. The constraint and objective vectors, matrices are set up with `numpy.ndarray` data types mostly and some basic NumPy functions. Thereafter `linprog()` from SciPy library is called alongwith integrality vectors via standard **HiGHS**(see [HiGHS documentation](https://ergo-code.github.io/HiGHS/dev/interfaces/cpp/library/)) solver to output the ILP value. Separately, the LP relaxation is also solved dropping the `integrality` vector to compute `result.ip_lp_ratio`.

~~~python
def generate_maximal_patterns(W, m, sizes):
~~~ 
This helper function generates maximal patterns for each bin corresponding to input instance `(W, m, sizes)`, still preserving the optimal objective of the normal Gilmore-Gomory Model additionally saving on number of varables/columns in the LP. I count the **number** of maximal patterns by **dynamic programming(DP)** in time $O(m W^2)$ and space $O(mW)$ and then the **list of all maximal patterns** in $O(m^2 P + m W P)$, where $P$ is the number of maximal patterns from the DP table.

`dp[i][r]` represents number of ways to coplete a maximal pattern using types $i,\dots,m - 1$, with remaining bin-capacity $r$. Then the recursion amounts to 
~~~python
dp[i][r] = sum([dp[i + 1][r - j * sizes[i]] for j in range(0, r // sizes[i] + 1)])
~~~

Now we compute the list of patterns using DP table with following idea: process items from $0$ to $m-1$. If say $[p_1, ..., p_z]$ is the (uncompleted) list of patterns at the $(i - 1)$ th stage, then for the ith stage see which $q$(= number of items) are feasible for each $p_k$ i.e. doesn't overshoot remaining capacity : `q * widths[i] <= remaining_capacity`, and at the same time corresponds( by extension ) to a maximal pattern via `dp[i + 1][remaining_capacity] > 0` , then append every such $q$ to $p_k$ to extend $p_k$ to (possibly)multiple patterns. Continue.

~~~python
def ggomory(W, m, sizes, time_limit):
~~~
This functions implements the Gilmore-Gomory model with maximal patterns, given `W`(bin capacity), `m`(number of items), `sizes`(the item size list) and `time_limit`(the time limit in seconds per instance). Similar in structure to the Kantorovich function.
~~~python
def arc_reduction1(W, out_adj, in_adj):
~~~
This helper function helps in reducing the number of arcs given the `W`(bin capacity), `out_adj`(the out-adjacencies in the arc-flow model graph), and `in_adj`(the in-adjacencies in the arc-flow model graph) according to the **Size ordering** reduction heuristic. Both `out_adj` and `in_adj` are maintained as a dictionary of python sets for $O(1)$ amortized access and deletion.
~~~python
def arc_reduction2(W, m, out_adj, in_adj, sizes, demands):
~~~
This helper function helps in reducing the number of arcs given the `W`(bin capacity),`m`(number of items), `out_adj` (the out-adjacencies in the arc-flow model graph), and `in_adj`(the in-adjacencies in the arc-flow model graph), `sizes`(the item sizes list) and `demands`(the item demands list) according to the **Demand-based reduction** reduction heuristic. Both `out_adj` and `in_adj` are maintained as a dictionary of python sets for $O(1)$ amortized access and deletion.

For each item size $i$, it collects all possible anchor points $P_i$, and for each $p \in P_i$ it collects all possible starting points $S_i^{p}$ for item-size $i$ arcs **under the reduction rule**. Thereafter it removes the arcs(of size $i$) which do not start at any point in $S_i = \cup_{p \in P_i}S_i^p$. **NOTE** : The **initial loss-arc removal** is performed implicitly during the arc-flow model graph construction and NOT implemented as another helper function.

~~~python
def arcflow(W, m, sizes, time_limit):
~~~
This functions implements the Arc-Flow model given `W`(bin capacity), `m`(number of items), `sizes`(the item size list) and `time_limit`(the time limit in seconds per instance).
It builds the graph by iterating over all vertices $0,\dots,W$ and at vertex $i$ adding all possible out-neighbors(at the same time, the correponding arcs in `in_adj`) corresponding to different sizes provided they don't exceed $W$.

Then it constructs the constraint matrices and contraint/objective vectors, placing arc $(i,j)$ at the `cumulative_degree_sum_list[i] + out_adj_list[i].index(j)`-th column, where `cumulative_degree_sum_list[i]` is the sum of outdegrees of nodes untill and including node $i-1$ and the flow/number-of-bins-used variable $z$ at the last column .

Thereafter, we check the input model and for each file-instance we run the model on the instance with the input time limit, and output 
```python
print(f"Instance {i + 1}/{num_instances}: "f"{result.message}", flush=True)
print("Construction + ILP solving:", result.total_time, "seconds") 
print("Only ILP solving:", result.ilp_time, "seconds")
```
for each case. If the ILP solve succeeds, we also output
```python
print(f"ILP optimum: {round(result.fun)}")
print("IP/LP:", result.ip_lp_ratio)
```

Finally for the whole folder/benchmark set, we output:
```python
print(f"Total instances solved to optimality: {instances_solved} out of {num_instances}")
print(f"Average total time: {avg_total_time / num_instances} seconds")
print(f"Average ILP time: {avg_ilp_time / num_instances} seconds")
if instances_solved == 0:
    print(f"Average bound strength: {None}")
else:
    print(f"Average bound strength: {avg_bound_strenth /      instances_solved}")
```
## Experimental evaluation
The report for the same is available [here](https://tilburgu-my.sharepoint.com/:b:/g/personal/a_gayen_tilburguniversity_edu/IQDYGT4gfoMlTKBCHKS13vOXAZ4N6zHK2tQ70qi3BsjtnPw?e=LV6wwb).


----
#### For feedback, corrections and  questions please send an email to `a.gayen@tilburguniversity.edu`. 
### Thank you for reading!



                                           





