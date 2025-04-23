Parameters
0        Model (0 = Set Partitioning; 1 = Set Covering) 
14400    Time Limit
2000000  Maximum number of columns
100000   Maximum number of columns for MIP
20.0     Initial alpha for the subgradient [20]
200      Maximum number of iterations without improvement for the subgradient [200]
0.00001  MipGap for the subgradient [0.00001]
0.25     PercK used to define how many columns generate at each iteration: k = PercK * n (where n is the number of items)  [0.25]
0.01     PercDelta: Delta = PercDelta * LB (Delta is the maximum reduced cost of the columns included in the core) 
0        Problem Type: 0=generate (using "Num_instances  Pro_Type  Num_items  Bin_Size  Constraint"; 1=read from "Num_instances  Filename"
Num_Istances
20
Num_instances  Pro_Type  Num_items  Bin_Size  Constraint
10                0          100       100         0 
10                0          200       100         0 
10                0          500       100         0 
10                0         1000       100         0 
10                0          100       500         0 
10                0          200       500         0 
10                0          500       500         0 
10                0         1000       500         0 
10                0          100      1000         0 
10                0          200      1000         0 
10                0          500      1000         0 
10                0         1000      1000         0 
10                0          100     10000         0 
10                0          200     10000         0 
10                0          500     10000         0 
10                0         1000     10000         0 
10                0          100    100000         0 
10                0          200    100000         0 
10                0          500    100000         0 
10                0         1000    100000         0 

