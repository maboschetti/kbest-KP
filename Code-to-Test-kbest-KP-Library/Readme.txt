++++++++++++  Code and executable  ++++++++++++ 

Folder "Code" includes the files:

kBest-KP-Library.cpp
kBest-KP-Library.h
Main-KP-Test.cpp
Utilities.cpp
Utilities.h


By running the code one should obtain an executable and call it KP-kbestLib-v3.exe 
This xecutable should be copied int both the "Tests" and "TestsLargeInstances" folder.


++++++++++++  Algorithm Name  ++++++++++++ 

The following use abbreviations for the algorithm types:

ALGORITHM KPFW_kbest_L is referred to as A4
ALGORITHM KPBW_kbest_L is referred to as A7 
ALGORITHM KPBW_kbest_LDCF(M) is referred to as A10

++++++++++++  Instances  ++++++++++++

The "Instances" folder  Includes the instance files (both normal size and large size)
Each file name reports: n_w_M_D_run.txt

n: number of items

w: maximum capacity

M: solution size

D: problem dimension

run: run identifier


++++++++++++  Tests on normal instances ++++++++++++  

The "Tests" folder containes all batch files needed to run the tests.

It must include the executable: KP-kbestLib-v3.exe

To start the testing process, run the file "tests.bat"

The file "tests.bat" will call: "testsA4.bat", "testsA7.bat", "testsA10.bat" 

Each of these will execute the corresponding algorithm (A4, A7, A10) for different runs (0, 1, 2) and various values of k.

Each line in the called batch files executes one test in the following format:
 
In those batch files, each line calls an algorithm (4, 7 or 10) to solve an instances and it is as follows:
  
executable	instance_name	value_of_k	random_id	id	algorithm_type	run_number
	
id: A unique identifier assigned to each instance, typically ordered by instance size. It provides a consistent reference to the same instance across all runs.

random_id: A run-specific identifier that randomly assigns an order to the instances for each run. This allows the instances to be processed in a different sequence across multiple runs, which is useful for analyzing performance variability or avoiding biases due to fixed execution order.


A "Results" folder is used to store the results on .cvs files, organized by algorithm type, k value, and number of the run.

The output file for the algorithm A4 shows these values:
-the name of the instance
-the algorithm type
-the id of the instance
-the id of the instance randomly generated
-the number of the run (0, 1, 2)
-n
-W
-k
-M
-D
-the best solution
-the kth solution (or the last solution found if k solutions do not exist)
-the generation computing time
-the total computing time

The output files for the algorithms A7 and A10 show these values:
-the name of the instance
-the algorithm type
-the id of the instance
-the id of the instance randomly generated
-the number of the run (0, 1, 2)
-n
-W
-k
-M
-D
-the best solution
-the kth solution (or the last solution found if k solutions do not exist)
-the average k
-the number of states used
-the % of states used
-the number of states generated
-the % of states used
-the average dimension of k
-the % dimension of k
-the generation computing time
-the total computing time


++++++++++++  Tests on large instances ++++++++++++  

The "TestsLargeInstances" folder contains bathc files on large intstnces sets.

It must include the executable: KP-kbestLib-v3.exe

To run the tests, execute "testsA10_large.bat". Other batch files function similarly to those used for normal instances.

Results will be stored in the "Results" folder as .csv files, organized by algorithm type, k value, and run number.

The output files for the algorithm A10 shows these values:
-the name of the instance
-the algorithm type
-the id of the instance
-the id of the instance randomly generated
-the number of the run (0, 1, 2)
-n
-W
-k
-M
-D
-the best solution
-the kth solution (or the last solution found if k solutions do not exist)
-the average k
-the number of states used
-the % of states used
-the number of states generated
-the % of states used
-the average dimension of k
-the % dimension of k
-the generation computing time
-the total computing time


